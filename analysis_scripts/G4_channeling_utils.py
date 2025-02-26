# Set of functions used during the analysis of channeling simulations with Geant4.

# Import the required libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import uproot
import math
from tqdm import tqdm
import scipy
from scipy.special import voigt_profile
from scipy.special import erf
import json

from G4_utils import * #my custom functions of general utility


#######################################################################################################
############ Custom functions useful for various Geant4 simulations related to channeling #############
def get_photons_on_detector(filename, Nevents, xlim_rad=(0, 1e10), \
                            apply_collimation=False, cut_angle=3.14, \
                            thetaC=0, beVerbose=True):
    """
    Function to open the root file obtained with the Geant4 FastSimChannelingRad app 
    and get the data scored in the "photon_spectrum" ntuple, namely the features of the 
    photons that impinge as the detector: ["E", "angle_x", "angle_y", eventID], with 
    E is expressed in MeV and the angles in rad. 
    The function returns the read root file and a dataframe with the selected photons.
    The angular selection can be applied specifying a cut angle with respect to a center thetaC.
    """
    import numpy as np
    import pandas as pd
    import uproot
    
    rf = uproot.open(filename)
    rf_content = [item.split(';')[0] for item in rf.keys()]
    print('rf_content:\n', rf_content, '\n')
    
    # get the simulated data
    if 'photon_spectrum' in rf_content:
        df_ph = rf['photon_spectrum'].arrays(library='pd')
    else:
        #Evalues = np.zeros(1)
        df_root = rf['scoring_ntuple'].arrays(library='pd')
        df_det_ph = df_root[(df_root.volume == "Detector") & (df_root.particle == "gamma")].copy()
        df_ph = df_det_ph.drop(["particleID", "parentID", "particle", "volume"], axis=1)
    Evalues = df_ph["E"].values #MeV
    print("number of photons scored:", Evalues.shape[0])
    
    # Take only the photons inside the collimator acceptance
    thetaX = df_ph["angle_x"].values #rad
    thetaY = df_ph["angle_y"].values #rad
    df_ph["theta"] = np.sqrt(thetaX**2 + thetaY**2)
    if apply_collimation:
        df_ph_sel = df_ph[np.abs(df_ph["theta"] - thetaC) < cut_angle]
    else:
        df_ph_sel = df_ph.copy()    
    df_ph_sel = df_ph_sel[(df_ph_sel.E >= xlim_rad[0]) & (df_ph_sel.E <= xlim_rad[1])]
    eventID_sel = list(df_ph_sel.eventID)

    # Print number of photons emitted
    if beVerbose:
        print("number of collimated photons:", len(df_ph[np.abs(df_ph["theta"] - thetaC) < cut_angle]))
        print("number of photons emitted within %.3f mrad with energy in [%.2f, %.2f] MeV: %d\n" % \
              (cut_angle*1e3, *xlim_rad, len(df_ph_sel.E)))
        print("\nnumber of photons emitted per particle: %.2f" % (len(df_ph)/Nevents))
        print("number of photons emitted per particle within %.3f mrad with energy in [%.2f, %.2f] MeV: %.4f\n" % \
              (cut_angle*1e3, *xlim_rad, len(df_ph_sel.E)/Nevents))

    return rf, df_ph_sel


def read_G4_BK_spectrum(filename, th=1.):
    """
    Function to read a text file with the photon spectrum emitted by an oriented crystal,
    analitically calculated in Geant4 through Baier-Katkov formula using virtual photons. 
    It returns [Esteps, spectrum, Nevents, Nspectra, Nbroken, Ne].
    """
    
    weights_list = []
    partial = []
    data_read = []
    broken = False
    Nspectra = 0
    Nbroken = 0

    # Read the file and create a dataframe with all the good (value < th) file entries.
    # NOTE: after last modifications, broken data are not present anymore.
    with open(filename, 'r') as f:
        header = f.readline()
        weights_list.append(int(header.split('\n')[0]))
        for line in f:
            if not line == '\n': 
                #E = float(line.split(' ')[0]) #MeV
                dNdE = float(line.split(' ')[1].split('\n')[0])
                #EdNdE = float(line.split(' ')[2].split('\n')[0])           
                partial.append(line)               
                if dNdE > th: #check using the threshold
                    broken = True    
            else:    
                if not broken:
                    data_read.append(partial)
                    Nspectra += 1
                else:
                    Nbroken += 1
                header = f.readline()
                if not header == "":
                    weights_list.append(int(header.split('\n')[0]))
                if broken:
                    weights_list.pop()
                broken = False
                partial = []
        if not broken:
            data_read.append(partial)  

    weights = np.array(weights_list, dtype='float')
    data_read = [item for sublist in data_read for item in sublist]
    Ndata = len(data_read)
    #print('weights_list:', weights_list)
    #print(data_read)
    #print("Nbroken:", Nbroken)
    #print("Nspectra:", Nspectra)

    data_x_read = np.array([data_read[i].split(' ')[0] for i in range(Ndata)], dtype='float')
    data_y_read = np.array([data_read[i].split(' ')[1].split('\n')[0] for i in range(Ndata)], dtype='float')
    data_z_read = np.array([data_read[i].split(' ')[2].split('\n')[0] for i in range(Ndata)], dtype='float')

    df_list = np.transpose(np.array([data_x_read, data_y_read, data_z_read]))
    df = pd.DataFrame(data=df_list, columns=['E', 'dNdE', 'EdNdE'])
    #print('\n', df.head(5), '\n')

    # build the spectrum (and the spectral intensity)
    Esteps = np.unique(df.E.values)
    Nbin = Esteps.shape[0]
    spectrum = np.zeros(Nbin)
    spectral_intensity = np.zeros(Nbin)

    my_list = []
    for E in list(Esteps):
        temp = df.loc[df['E'] == E].values
        my_list.append(temp)    
    my_array = np.array(my_list, dtype='float')

    for i in range(Nbin):
        my_array2 = my_array[i]
        spectrum[i] = np.average(my_array2[:,1], weights=weights)
        spectral_intensity[i] = np.average(my_array2[:,2], weights=weights)   

    spectrum[-1] = 0
    
    Nevents = int(np.sum(weights))
    Ne = round(np.sum(spectrum), 2)
    print('Number of good events:', Nevents)
    print('Number of photons emitted per primary (single-photon approx):', Ne)
    
    return Esteps, spectrum, Nevents, Nspectra, Nbroken, Ne


def get_deflection_angles(df_out, df_in, pot_good_events, 
                          res_thetaX_in=0, res_thetaY_in=0, \
                          res_DthetaX=0, res_DthetaY=0, \
                          ang_cut=np.pi*1e6, cut_center=[0., 0.]):
    """
    Calculate the deflection for events whose transverse
    deviation angle is within a given angular cut.
    Experimental resolution can be taken into account.
    Input angles must be given in [urad]. However, the angles
    in the passed dataframes (df_out, df_in) are expressed in [rad].
    """
    DthetaX = [] #urad
    DthetaY = [] #urad
    thetaX_in = [] #rad
    thetaY_in = [] #rad
    good_events = []
    for i in tqdm(pot_good_events):
        df_out_i = df_out[df_out["eventID"] == i].copy()
        df_in_i = df_in[df_in["eventID"] == i].copy()  
        thX_in = np.random.normal(df_in_i['angle_x'].values[0]*1e6, res_thetaX_in, 1)[0] #urad
        thY_in = np.random.normal(df_in_i['angle_y'].values[0]*1e6, res_thetaY_in, 1)[0] #urad
        if np.sqrt((thX_in-cut_center[0])**2 + (thY_in-cut_center[1])**2) < ang_cut: #urad
            thetaX_in.append(thX_in)
            thetaY_in.append(thY_in)
            DthetaX_mean = (df_out_i['angle_x'].values[0] - df_in_i['angle_x'].values[0])*1e6 #urad
            DthetaY_mean = (df_out_i['angle_y'].values[0] - df_in_i['angle_y'].values[0])*1e6 #urad  
            DthetaX.append(np.random.normal(DthetaX_mean, res_DthetaX, 1)[0])
            DthetaY.append(np.random.normal(DthetaY_mean, res_DthetaY, 1)[0])
            good_events.append(i)
        del df_out_i, df_in_i
    return DthetaX, DthetaY, thetaX_in, thetaY_in, good_events


def filter_photon_emission_with_defl_cut(DthetaX, DthetaY, good_events, df_ph, \
                                         applySel=False, selType=1, \
                                         DthXmin=-np.pi*1e6, DthXmax=np.pi*1e6, \
                                         DthYmin=-np.pi*1e6, DthYmax=np.pi*1e6, \
                                         DthCx=0, DthCy=0, \
                                         DthR=2*np.pi*1e6, \
                                         DthCx2=0, DthCy2=0, \
                                         DthR2=2*np.pi*1e6, \
                                         DthCx3=0, DthCy3=0, \
                                         DthR3=2*np.pi*1e6, \
                                        ):
    """
    Filter the photon emission from a crystal based on the cuts on deflection. A squared (selType=0)
    or circular (selType=1) cut with parameters passed as arguments (in urad) can be applied.
    I updated the function to consider also 2 or 3 circular selctions, as in the BentGe<110> article. 
    DthetaX, DthetaY, good_events are lists or numpy arrays. good_events lists some events that were 
    previously selected according to some criteria (i.e. previous angular/spatial cuts).
    df_ph is the dataframe of emitted photons. It contains ['eventID', 'Ekin', 'angle_x', 'angle_y'].
    The function returns both the total energy lost by radiation per each selected event (Erad_sel),
    and the individual photon emission features (Eph_sel[GeV], thetaX_ph_sel[mrad], thetaY_ph_sel[mrad])
    as numpy arrays.
    """
    Eph_dict_sel = {} #it will contain the energies of photons emitted at each event
    thetaX_ph_dict_sel = {} #it will contain the angle_x of photons emitted at each event
    thetaY_ph_dict_sel = {} #it will contain the angle_y of photons emitted at each event
    sel_events = []
    if applySel:
        print("selType:", selType)
        if selType == 0:
            for i in range(len(good_events)):
                if (DthXmin <= DthetaX[i] and DthetaX[i] <= DthXmax) and \
                   (DthYmin <= DthetaY[i] and DthetaY[i] <= DthYmax):
                    sel_events.append(good_events[i])
        elif selType == 1:
            for i in range(len(good_events)):
                if ( (DthetaX[i] - DthCx)**2 + (DthetaY[i] - DthCy)**2 <= DthR**2 ):
                    sel_events.append(good_events[i])
        elif selType == 2:
            for i in range(len(good_events)):
                if ( ((DthetaX[i] - DthCx2)**2 + (DthetaY[i] - DthCy2)**2 <= DthR2**2) or \
                     ((DthetaX[i] - DthCx3)**2 + (DthetaY[i] - DthCy3)**2 <= DthR3**2) ):
                    sel_events.append(good_events[i])
        else:
             for i in range(len(good_events)):
                if not ( ((DthetaX[i] - DthCx)**2 + (DthetaY[i] - DthCy)**2 <= DthR**2) or \
                         ((DthetaX[i] - DthCx2)**2 + (DthetaY[i] - DthCy2)**2 <= DthR2**2) or \
                         ((DthetaX[i] - DthCx3)**2 + (DthetaY[i] - DthCy3)**2 <= DthR3**2) ):
                    sel_events.append(good_events[i]) 
    else:
        sel_events = good_events
    for i in range(len(sel_events)):
        df_i = df_ph[df_ph.eventID == sel_events[i]]
        Eph_dict_sel[sel_events[i]] = list(df_i["Ekin"].values*0.001) #MeV -> GeV
        thetaX_ph_dict_sel[sel_events[i]] = list(df_i["angle_x"].values*1e3) #rad -> mrad
        thetaY_ph_dict_sel[sel_events[i]] = list(df_i["angle_y"].values*1e3) #rad -> mrad
    Erad_sel = np.array([sum(Eph_dict_sel[event]) for event in Eph_dict_sel.keys()]) #total energy lost by radiation per (selected) event [GeV]
    print("number of selected events:", len(Eph_dict_sel.keys()))
    Eph_sel = [] #it will contain the individual photon spectrum [GeV]
    thetaX_ph_sel = [] #it will contain the individual photon angle_x [mrad]
    thetaY_ph_sel = [] #it will contain the individual photon angle_y [mrad]
    temp1 = []
    temp2 = []
    temp3 = []
    for event in Eph_dict_sel.keys():
        temp1 = Eph_dict_sel[event]
        Eph_sel.append(temp1)
        temp1 = []
        temp2 = thetaX_ph_dict_sel[event]
        thetaX_ph_sel.append(temp2)
        temp2 = []
        temp3 = thetaY_ph_dict_sel[event]
        thetaY_ph_sel.append(temp3)
        temp3 = []
    Eph_sel = np.array(list_flatten(Eph_sel))
    thetaX_ph_sel = np.array(list_flatten(thetaX_ph_sel))
    thetaY_ph_sel = np.array(list_flatten(thetaY_ph_sel))
    print("number of selected photons:", len(Eph_sel) , '\n')
    return Erad_sel, Eph_sel, thetaX_ph_sel, thetaY_ph_sel
#######################################################################################################

#######################################################################################################
#### Functions implemented by R. Negrello to fit theta_out curve considering various contributions ####
path_file_json = "parameters.json"
if os.path.exists(path_file_json):
    with open(path_file_json, 'r') as file:
        parameters = json.load(file)

    def f_dech_exp(theta_X, theta_dech, A_dech):
        sigma_VR = parameters["sigma_VR"]  
        theta_VR = parameters["mu_VR"]
        theta_Ch = parameters["mu_ch"]
        term1 = A_dech / (2 * theta_dech)
        term2_numerator = (sigma_VR**2) / (2 * (theta_dech**2)) + (theta_Ch - theta_X) / theta_dech
        term2 = np.exp(term2_numerator)
        term3_numerator1 = theta_VR - theta_X + (sigma_VR**2) / theta_dech
        term3_denominator = np.sqrt(2) * sigma_VR
        term3_erf1 = erf(term3_numerator1 / term3_denominator)
        term3_numerator2 = theta_Ch - theta_X + (sigma_VR**2) / theta_dech
        term3_erf2 = erf(term3_numerator2 / term3_denominator)
        result = term1 * term2 * (term3_erf1 - term3_erf2)
        return result

    def f_dech_sim(theta_X, theta_dech, A_dech):
        sigma_VR = parameters["sigma_VR_sim"]  
        theta_VR = parameters["mu_VR_sim"]
        theta_Ch = parameters["mu_ch_sim"]
        term1 = A_dech / (2 * theta_dech)
        term2_numerator = (sigma_VR**2) / (2 * (theta_dech**2)) + (theta_Ch - theta_X) / theta_dech
        term2 = np.exp(term2_numerator)
        term3_numerator1 = theta_VR - theta_X + (sigma_VR**2) / theta_dech
        term3_denominator = np.sqrt(2) * sigma_VR
        term3_erf1 = erf(term3_numerator1 / term3_denominator)
        term3_numerator2 = theta_Ch - theta_X + (sigma_VR**2) / theta_dech
        term3_erf2 = erf(term3_numerator2 / term3_denominator)
        result = term1 * term2 * (term3_erf1 - term3_erf2)
        return result
    
def proj(img, axis):
    axisN = 1 if axis=="x" else 0
    return np.sum(img.astype("float"), axisN)

def gauss(x, A, mu, sigma):
    return A/(sigma*np.sqrt(2*np.pi))*np.exp(-(x-mu)**2/(2*sigma**2))

def f_VR(x, A, mu, sigma, r):
    return A/(sigma*np.sqrt(2*np.pi))*np.exp(-(x-mu)**2/(2*sigma**2))+ \
          (1-A)/(r*sigma*np.sqrt(2*np.pi))*np.exp(-(x-mu)**2/(2*r**2*sigma**2))

def voigt(x, A, mu, sigma, gamma):
    return A/(sigma*np.sqrt(2*np.pi)) * voigt_profile(x - mu, sigma, gamma)    
#######################################################################################################

#######################################################################################################
################################ Functions for particle tagging #######################################
def my_dataframe_split(df, events, len_traj, Nparts_desired=500):
    """
    Custom method (by gpaterno) to split a dataframe in chunks/parts on the basis of
    a column (eventID), so as to avoid of having partial trajectories in each part.
    It returns a list of dataframes.
    """
    # preliminary operations
    Nevents = len(events)
    Nentries = len(df)
    max_traj_len = max(len_traj.values())
    Nparts_approx = min([int(Nentries / max_traj_len), Nparts_desired])
    NtrajXpart = int(Nevents / Nparts_approx)
    print("Nevents:", Nevents, ", Nentries:", Nentries)
    print("max trajectory length:", max_traj_len)
    print("NtrajXpart:", NtrajXpart)
    # calculation loop
    df_parts = []
    temp = []
    k = 0
    ke = 0 #number of processed events
    kp = 0 #number of created parts
    kr = 0 #number of residual events
    temp_r = []
    Nentries_processed = 0
    #for i in tqdm(events):
    #    df_i = df[df["eventID"] == i].copy() #variables of this event
    df_grouped = df.groupby('eventID')
    for i, df_i in tqdm(df_grouped):
        if i in events:
            temp.append(df_i)
            ke += 1
            k += 1            
            if k >= NtrajXpart:
                df_temp = pd.concat(temp, ignore_index=True)
                df_parts.append(df_temp)
                temp = []
                k = 0
                Nentries_processed += len(df_temp) 
                kp += 1           
            Nevents_saved = kp*NtrajXpart
            if (ke > Nevents_saved) and (Nevents - Nevents_saved < NtrajXpart):
                #print("ke", ke)
                temp_r.append(df_i)
                kr += 1
    df_temp = pd.concat(temp_r, ignore_index=True)
    df_parts.append(df_temp)
    Nentries_processed += len(df_temp)    
    # print
    #print("kp", kp) 
    #print("kr", kr) 
    print("df dataframe actually split into %d parts:" % len(df_parts))
    print("Nevents_processed:", ke)
    print("Nentries_processed:", Nentries_processed)
    if (Nentries - Nentries_processed == 0):
        print("All of those available!")
    else:
        print("WARNING: We miss %d entries" % (Nentries - Nentries_processed))
    # return a list of partial dataframes
    return df_parts


# Functions for parallel calculation of particle state (by R. Negrello with my mod)
def calculate_E_K(df_part, BeamEnergy):
    E_K_partial = []
    for _, row in df_part.iterrows():
        E_K = BeamEnergy/2 * row['tx']**2
        E_K_partial.append(E_K)
    return E_K_partial


def calculate_E_T(df_part, x, BeamEnergy, U_eff):
    E_T_partial = []
    for _, row in df_part.iterrows():
        i = np.abs(x - row['x']*1e7).argmin()
        E_T = BeamEnergy/2 * row['tx']**2 + U_eff[i]
        E_T_partial.append(E_T)
    return E_T_partial


def determine_states(df_part, x, BeamEnergy, U_eff, hkl, d_pl, particleCharge):
    states_partial = []
    if particleCharge > 0:
        toll = 150
        for _, row in df_part.iterrows():
            i = np.abs(x - row['x']*1e7).argmin()
            E_T = BeamEnergy/2 * row['tx']**2 + U_eff[i]
            if hkl == (1,1,1): #gpaterno
                if ((row['x']*1e7 >= 0) and (row['x']*1e7 <= d_pl[0]/2)):
                    if E_T <= max(U_eff[np.abs(-(d_pl[0]/2)-x).argmin()-toll:np.abs(-(d_pl[0]/2)-x).argmin()+toll]):
                        state = 'c'
                    else:
                        state = 'o'
                elif ((row['x']*1e7 >= d_pl[0]/2) and (row['x']*1e7 <= d_pl[0]/2+d_pl[1])):
                    if E_T <= max(U_eff[np.abs(d_pl[0]/2-x).argmin()-toll:np.abs(d_pl[0]/2-x).argmin()+toll]):
                        state = 'c'
                    else:
                        state = 'o'
                else:
                    if E_T <= max(U_eff[np.abs((d_pl[0]/2+d_pl[1])-x).argmin()-toll:np.abs((d_pl[0]/2+d_pl[1])-x).argmin()+toll]):
                        state = 'c'
                    else:
                        state = 'o'
            else: #(1,1,0) by gpaterno, to check!
                if ((row['x']*1e7 >= 0) and (row['x']*1e7 <= d_pl[0]/2)):
                    if E_T <= (U_eff[np.abs(-d_pl[0]/2-x).argmin()]):
                        state = 'c'
                    else:
                        state = 'o' 
                else:
                    if E_T <= (U_eff[np.abs(d_pl[0]/2-x).argmin()]):
                        state = 'c'
                    else:
                        state = 'o'                   
            states_partial.append(state)
    else:
        toll = 50
        for _, row in df_part.iterrows():
            i = np.abs(x - row['x']*1e7).argmin()
            E_T = BeamEnergy/2 * row['tx']**2 + U_eff[i]
            if E_T <= max(U_eff[np.abs(0 - x).argmin()-toll:np.abs(0 - x).argmin()+toll]):
                state = 'c'
            else:
                state = 'o'
            states_partial.append(state)
    return states_partial


def parallelize_dataframe_parts(df_parts, func, *args, Nproc_des=-1):
    import multiprocessing
    if Nproc_des == -1: #gpaterno
        num_processes = multiprocessing.cpu_count()
        print("num_processes:", num_processes)
    else:
        num_processes = Nproc_des
    pool = multiprocessing.Pool(num_processes)
    results = pool.starmap(func, [(df_part, *args) for df_part in df_parts])
    pool.close()
    pool.join()
    return results


def classify_particles(df_part):
    """
    Function to classify particle trajectories with or without parallel calculation 
    (by R. Negrello). It returns a dataframe with two columns [eventID', 'category']. 
    To avoid future disambiguations, I changed the second column name from 'state' to
    'category' and the returned dataframe name from 'df_state' to 'df_category'.
    In order to correctly use this function, include here (or before calling):
    from collections import Counter
    tag_counters = Counter()
    """
    
    from collections import Counter
    tag_counters = Counter()
    
    df_category = pd.DataFrame(columns=['eventID', 'category'])

    for event_id, group in tqdm(df_part.groupby(['eventID'])):
        
        states = group['state'].tolist()
        final_state = 'Unknown'
        
        # Counts state occurrence
        state_counts = Counter(states)
        
        # Channeling
        if state_counts['c']==len(states) and states[0]=="c" and states[-1]=="c": 
            final_state = 'ch'
            
        # Overbarrier
        elif state_counts['o']==len(states):
                final_state = 'ob'        
        else:

            if "o" in states and "c" in states:
                ob_index = states.index('o')
                ch_index = states.index('c')  

            # Captured Overbarrier
                if states[0]=="o": # and all(state == 'o' for state in states[:ch_index]):
                    if all(state=='c' for state in states[ch_index:]):
                        final_state = 'co+ch'
                    elif all(state in ['c', 'o'] for state in states[ch_index:]):
                        final_state = 'co+de'
                
            # Dechanneling 
                elif (states[0]=="c") and (states[-1]=="o") \
                                      and all(state=='c' for state in states[:ob_index]) \
                                      and all(state=='o' for state in states[ob_index:]):
                        final_state = 'dech'
                
                # Rechanneling
                else:
                    ch_indexes = [i for i, state in enumerate(states) if state=='c']
                    num_rech_transitions = 0
                    
                    # check transitions from "ch" to "non-ch" states and then again to "ch"
                    for i in range(len(ch_indexes) - 1):
                        start_index = ch_indexes[i]
                        end_index = ch_indexes[i+1]
                        if any(state!='c' for state in states[start_index+1:end_index]) and states[end_index]=='c':
                            num_rech_transitions += 1
                            
                    if num_rech_transitions == 1:
                        if states[-1]=="c":
                            final_state = 'rech1'
                        else:
                            final_state = 'rech1+dech'
                    elif num_rech_transitions == 2:
                        if states[-1]=="c":
                            final_state = 'rech2'
                        else:
                            final_state = 'rech2+dech'
                    elif num_rech_transitions == 3:
                        if states[-1]=="c":
                            final_state = 'rech3'
                        else:
                            final_state = 'rech3+dech'
                    elif num_rech_transitions == 4 :
                        if states[-1]=="c":
                            final_state = 'rech4'
                        else:
                            final_state = 'rech4+dech'
                    elif num_rech_transitions == 5 :
                        if states[-1]=="c":
                            final_state = 'rech5'
                        else:
                            final_state = 'rech5+dech'
                    elif num_rech_transitions == 6 :
                        if states[-1]=="c":
                            final_state = 'rech6'
                        else:
                            final_state = 'rech6+dech'
                    elif num_rech_transitions == 7:
                        if states[-1]=="c":
                            final_state = 'rech7'
                        else:
                            final_state = 'rech7+dech'
                    elif num_rech_transitions == 8:
                        if states[-1]=="c":
                            final_state = 'rech8'
                        else:
                            final_state = 'rech8+dech'
                    elif num_rech_transitions == 9 :
                        if states[-1]=="c":
                            final_state = 'rech9'
                        else:
                            final_state = 'rech9+dech'
                    elif num_rech_transitions == 10 :
                        if states[-1]=="c":
                            final_state = 'rech10'
                        else:
                            final_state = 'rech10+dech'
                    elif num_rech_transitions == 11 :
                        if states[-1]=="c":
                            final_state = 'rech11'
                        else:
                            final_state = 'rech11+dech'
                    elif num_rech_transitions == 12:
                        if states[-1]=="c":
                            final_state = 'rech12'
                        else:
                            final_state = 'rech12+dech'
                    elif num_rech_transitions == 13:
                        if states[-1]=="c":
                            final_state = 'rech13'
                        else:
                            final_state = 'rech13+dech'
                    elif num_rech_transitions == 14 :
                        if states[-1]=="c":
                            final_state = 'rech14'
                        else:
                            final_state = 'rech14+dech'
                    elif num_rech_transitions == 15 :
                        if states[-1]=="c":
                            final_state = 'rech15'
                        else:
                            final_state = 'rech15+dech'
                    elif num_rech_transitions == 16 :
                        if states[-1]=="c":
                            final_state = 'rech16'
                        else:
                            final_state = 'rech16+dech'
                    elif num_rech_transitions == 17:
                        if states[-1]=="c":
                            final_state = 'rech17'
                        else:
                            final_state = 'rech17+dech'
                    elif num_rech_transitions == 18:
                        if states[-1]=="c":
                            final_state = 'rech18'
                        else:
                            final_state = 'rech18+dech'
                    elif num_rech_transitions == 19 :
                        if states[-1]=="c":
                            final_state = 'rech19'
                        else:
                            final_state = 'rech19+dech'
                    elif num_rech_transitions == 20 :
                        if states[-1]=="c":
                            final_state = 'rech20'
                        else:
                            final_state = 'rech20+dech'

        df_category = pd.concat([df_category, \
                                 pd.DataFrame({'eventID': event_id, 'category': [final_state]})], \
                                              ignore_index=True)
    return df_category


def correct_trajectory_of_radiating_particles(df):
    """
    Function, initially developed by R. Negrello and generalized by gpaterno, 
    to correct the trajectory of the channeled particles that radiate. Indeed, 
    in this case, the coordinates of the trajectory invert one or more times. 
    The selection of the radiating particles is carried out a-priori outside 
    of the function. Indeed, a filtered version of the dataframe read from the
    "tag" txt file is passed as input. The same dataframe, but with 
    variables properly corrected, is returned. 
    """   

    df_corr = pd.DataFrame(columns=df.columns)
    df_grouped = df.groupby('eventID')
    for i, df_i in tqdm(df_grouped):
        x_val = df_i['x'].values
        tx_val = df_i['tx'].values
        z_val = df_i['z'].values
        xx_val = df_i['xx'].values
        
        inds_inv_coord = []
        for i in range(len(z_val) - 1):
            if (z_val[i+1] < z_val[i]):
                inds_inv_coord.append(i+1)
        #print("inds_inv_coord of event %d:" % i, inds_inv_coord)

        emission_idxs = []
        for ii in inds_inv_coord:
            emission_idxs.append(np.abs(z_val[ii]-z_val[:ii]).argmin())
        couples = np.stack([inds_inv_coord, emission_idxs],axis=1)

        if len(inds_inv_coord) == 1:
            new_x_val = np.concatenate([x_val[:couples[0][1]], x_val[couples[0][0]:]])
            new_tx_val = np.concatenate([tx_val[:couples[0][1]], tx_val[couples[0][0]:]])
            new_z_val = np.concatenate([z_val[:couples[0][1]], z_val[couples[0][0]:]])

        if len(inds_inv_coord) == 2:
            new_x_val = np.concatenate([x_val[:couples[0][1]], x_val[couples[0][0]:couples[1][0]]])
            new_tx_val = np.concatenate([tx_val[:couples[0][1]], tx_val[couples[0][0]:couples[1][0]]])
            new_z_val = np.concatenate([z_val[:couples[0][1]], z_val[couples[0][0]:couples[1][0]]])

        if len(inds_inv_coord) == 3:
            new_x_val = np.concatenate([x_val[:couples[0][1]], x_val[couples[0][0]:couples[1][1]], x_val[couples[1][1]:couples[2][1]]])
            new_tx_val = np.concatenate([tx_val[:couples[0][1]], tx_val[couples[0][0]:couples[1][1]], tx_val[couples[1][1]:couples[2][1]]])
            new_z_val = np.concatenate([z_val[:couples[0][1]], z_val[couples[0][0]:couples[1][1]], z_val[couples[1][1]:couples[2][1]]])

        if len(x_val) > len(new_x_val):
            df_i['x'] = np.pad(new_x_val, (0, len(x_val)-len(new_x_val)), 'constant')
            df_i['tx'] = np.pad(new_tx_val, (0, len(tx_val)-len(new_tx_val)), 'constant')
            df_i['z'] = np.pad(new_z_val, (0, len(z_val)-len(new_z_val)), 'constant')
            meanX = np.mean(new_x_val)
            xarray = np.array(new_x_val)
            xarray = np.where(xarray < meanX, xarray + meanX, xarray - meanX)
            xarray -= meanX
            df_i['xx'] = np.pad(xarray, (0, len(xx_val)-len(xarray)), 'constant')

            df_corr = df_i.copy() if df_corr.empty else pd.concat([df_corr, df_i], ignore_index=True)
        else:
            print("particle %d discarded" % i)
    
    print("correction done!")
    return df_corr
#######################################################################################################

#######################################################################################################
################## Merge the results of many simulations with few (even 1) events #####################
def merge_ntuples_from_files_into_df(data_path, ntuple_name):
    """
    This function is quite general, it merges ntuples in python.
    There are no constranints on the filenames. It returns a dataframe.
    """
    import uproot
    dataframes = []
    for filename in os.listdir(data_path):
        print('opening', filename)
        file_path = os.path.join(data_path, filename)
        if os.path.isfile(file_path):
            rf = uproot.open(file_path)
            rf_content = [item.split(';')[0] for item in rf.keys()]
            if ntuple_name in rf_content:
                df = rf[ntuple_name].arrays(library='pd')
            dataframes.append(df)
    df_merged = pd.concat(dataframes, ignore_index=True)
    print('ntuples merged into a single dataframe!')
    return df_merged


def merge_files_into_dfs(data_path, beVerbose=False, save_result=False):
    """
    Function to merge both the root and the txt files of a set of
    Geant4 simulations of particle interaction in Oriented Crystals.
    The name of the corresponding root and txt files MUST be the same.
    It returns two dataframes, but it can also create merged files.
    """
    import uproot
    dataframes_defl_list = []
    dataframes_rad_list = []
    dataframes_txt_list = []
    events_read = 0
    nfiles = 0
    #IDs = [int(filename.split('.')[0]) for filename in os.listdir(data_path)]
    for filename in os.listdir(data_path):
        if beVerbose:
            print('opening', filename, '...')    
        # open root file
        file = os.path.join(data_path, filename)
        rf = uproot.open(file)
        rf_content = [item.split(';')[0] for item in rf.keys()]      
        # deflection ntuple
        df_defl = rf['scoring_ntuple'].arrays(library='pd')          
        # radiation ntuple
        df_rad = rf['photon_spectrum'].arrays(library='pd')
        df_rad = df_rad.rename(columns={"E": "Ekin"})            
        # list of events
        event_list = list(df_defl.eventID.unique())
        n_events = len(event_list)
        # correct the event number
        for event in event_list:
            df_defl.loc[df_defl.eventID == event, 'eventID'] += events_read
            dataframes_defl_list.append(df_defl)
            if not df_rad.empty:
                df_rad.loc[df_rad.eventID == event, 'eventID'] += events_read
                dataframes_rad_list.append(df_rad)           
        # open txt file
        volume = []
        eventID = []
        trackID = [] 
        x = []
        tx = []
        z = []
        xx = []
        yy = []
        file = os.path.join(data_path.replace('rad', 'tag'), filename.replace('root', 'txt'))     
        with open(file, 'r') as txt_file:
            for line in txt_file:
                values = line.strip().split()
                eventID.append(int(values[0]) + events_read)
                trackID.append(int(values[1]))
                x.append(float(values[2]))
                tx.append(float(values[3]))
                z.append(float(values[4]))
                xx.append(float(values[5]))
                #yy.append(float(values[6]))  
        df = pd.DataFrame({
            "volume": "Crystal",
            "eventID": eventID,
            "trackID": trackID,
            "x": x,
            "tx": tx,
            "z": z,
            "xx": xx
            #"yy": yy
        })
        dataframes_txt_list.append(df)
        # increase te number of files read
        events_read += n_events
        nfiles += 1
    # merge the dataframes
    df_defl_merged = pd.concat(dataframes_defl_list, ignore_index=True)
    df_defl_merged.particle = ["e-"]*len(df_defl_merged) #correction #<--- it may not be required!
    df_rad_merged = pd.concat(dataframes_rad_list, ignore_index=True)   
    df_txt_merged = pd.concat(dataframes_txt_list, ignore_index=True)
    # save merged daftrame to files
    if save_result:
        # export merged defl and rad dataframes to root files
        rf_merged = uproot.recreate(radfilespath + '../merged_root_file.root')
        tree_defl = rf_merged.mktree("scoring_ntuple", {
                                                 "eventID": np.int32, "volume": np.int32, \
                                                 "x": np.float64, "y": np.float64, \
                                                 "angle_x": np.float64, "angle_y": np.float64, \
                                                 "Ekin": np.float64, "particle": np.int32, \
                                                 "particleID": np.int32, "parentID": np.int32
                                                })
        df_defl_merged = df_defl_merged.replace({"Crystal": 0, "Detector": 1}) #without this an error occurs 
        df_defl_merged = df_defl_merged.replace({"e-": -1, "e+": 1, None: 0})  #without this an error occurs
        tree_defl.extend({
                          "eventID": df_defl_merged.eventID.values, "volume": df_defl_merged.volume.values, \
                          "x": df_defl_merged.x.values, "y": df_defl_merged.y.values, \
                          "angle_x": df_defl_merged.angle_x.values, "angle_y": df_defl_merged.angle_y.values, 
                          "Ekin": df_defl_merged.Ekin.values, "particle": df_defl_merged.particle.values, \
                          "particleID": df_defl_merged.particleID.values, "parentID": df_defl_merged.parentID.values
                         })
        tree_rad = rf_merged.mktree("photon_spectrum", {
                                                 "Ekin": np.float64, \
                                                 "angle_x": np.float64, "angle_y": np.float64, \
                                                 "eventID": np.int32,
                                                })    
        tree_rad.extend({
                         "Ekin": df_rad_merged.Ekin.values, \
                         "angle_x": df_rad_merged.angle_x.values, "angle_y": df_rad_merged.angle_y.values, \
                         "eventID": df_rad_merged.eventID.values
                        })
        # export merged txt dataframe to text file
        with open(radfilespath + '../merged_txt_file.txt', 'a') as f:
            df_string = df_txt_merged.to_string(header=False, index=False)
            f.write(df_string)
    # return merged dataframes
    print('%d files merged!\n' % (nfiles))
    return df_defl_merged, df_rad_merged, df_txt_merged
#######################################################################################################
