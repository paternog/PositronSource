#######################################################################################################
####### Set of functions used during the analysis of channeling simulations with Geant4 ###############
####### Author: Gianfranco Paternò (paterno@fe.infn.it), last update: 19/12/2025 ######################
#######################################################################################################

# Import the required libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import uproot
import math
from tqdm.notebook import tqdm
import scipy
from scipy.special import voigt_profile
from scipy.special import erf
import json
import time

from G4_utils import * #my custom functions of general utility


#######################################################################################################
############ Custom functions useful for various Geant4 simulations related to channeling #############
#######################################################################################################
def get_photons_on_detector(filename, Nevents, Elim=(0, 1e10), \
                            apply_collimation=False, coll_type='ellipse', tilt=45, \
                            cut_angle=(3.14, 3.14), thetaC=(0, 0), beVerbose=True):
    """
    Function to open the root file obtained with the Geant4 FastSimChannelingRad app 
    and get the data scored in the ntuples so as to select the main features of the photons 
    that impinge as the detector, namely: ["E", "angle_x", "angle_y", "eventID"], 
    with E is expressed in MeV and the angles in rad. 
    The function returns the read root file and a dataframe with the selected photons.
    The angular selection can be applied specifying cut angles with respect to a center thetaC.
    The collimator can be elliptical (circular) or rectangular. The selection can be performed
    using coll_type variable. Tilt and thetaC define the tilt of an elliptical collimator and
    the cut center, respectively. cut_angle define the abosute cut angles in x and y direction
    (for a circular collimator, the two elements must have the same value).
    Last update: 31/10/2025.
    """
    
    import numpy as np
    import pandas as pd
    import uproot
    
    rf = uproot.open(filename)
    rf_content = [item.split(';')[0] for item in rf.keys()]
    print('rf_content:\n', rf_content, '\n')
    
    ph_features = ["E", "angle_x", "angle_y", "eventID"]
       
    # Get the simulated data
    if 'photon_spectrum' in rf_content:
        df_ph = rf['photon_spectrum'].arrays(library='pd')
    elif 'detector_photons' in rf_content:
        df_root = rf['detector_photons'].arrays(library='pd')
        df_ph = [ph_features]
    else:
        df_root = rf['scoring_ntuple'].arrays(library='pd')
        df_det_ph = df_root[(df_root.volume == "Detector") & (df_root.particle == "gamma")].copy()
        df_ph = df_det_ph[ph_features]
    Evalues = df_ph["E"].values #MeV
    print("number of photons scored:", Evalues.shape[0])
       
    # Take only the photons inside the collimator acceptance
    if apply_collimation:
        if coll_type == 'ellipse': #or circle if cut_angle[0]=cut_angle[1]
            df_ph_sel = df_ph[np.sqrt((df_ph["angle_x"] - thetaC[0])**2/cut_angle[0]**2 + \
                                      (df_ph["angle_y"] - thetaC[1])**2/cut_angle[1]**2) < 1]
            #_, _, mask = elliptical_selection(df_ph["angle_x"], df_ph["angle_y"], thetaC, *cut_angle*2, tilt) #defined in G4_utils.py
            #df_ph_sel = df_ph[mask] #It's an equivalent method!
        else:
            mask1 = np.abs(df_ph["angle_x"] - thetaC[0]) < cut_angle[0]
            mask2 = np.abs(df_ph["angle_y"] - thetaC[1]) < cut_angle[1]
            df_ph_sel = df_ph[mask1 & mask2]
    else:
        df_ph_sel = df_ph.copy()    
    df_ph_sel_E = df_ph_sel[(df_ph_sel.E >= Elim[0]) & (df_ph_sel.E <= Elim[1])].copy()

    # Print number of photons emitted
    if beVerbose:
        print("number of collimated photons: %d" % len(df_ph_sel))
        print("number of collimated photons with energy in [%.2f, %.2f] MeV: %d" % \
              (*Elim, len(df_ph_sel_E.E)))
        print("number of photons emitted per particle: %.2f" % (len(df_ph)/Nevents))
        print("number of collimated photons emitted per particle with energy in [%.2f, %.2f] MeV: %.4f\n" % \
              (*Elim, len(df_ph_sel_E.E)/Nevents))
    
    # Return
    return rf, df_ph_sel, df_ph_sel_E


def read_G4_BK_spectrum(filename, file_format='new', th=1.):
    """
    Function to read a text file with the photon spectrum emitted by an oriented crystal
    analitically calculated in Geant4 through Baier-Katkov formula using sampling photons. 
    NOTE: the code is mainly useful to read a file contaning the results provided by 
    each thread separated by the number of events of the thread (file_format='old').
    Also, the code was developed when there was the possibility that the results of some 
    thread were bugged. For this reason there are many check, that actually are not needed anymore.
    If file_format='old' the code returns [Esteps, spectrum, Nevents, Nspectra, Nbroken, Ne].
    The new file format is much more straightforward to read, since it ha the format of a 
    text file with 2 columns: E_photon [MeV]  dW_rad/dE_photon [MeV^-1]. Threfore, it can
    be read simply with: np.loadtxt(filename, delimiter=' ', skiprows=1, unpack=True).
    Last update: 29/11/2025.
    """

    if file_format == 'old':
        weights_list = []
        partial = []
        data_read = []
        broken = False
        Nspectra = 0
        Nbroken = 0
    
        # Read the file and create a dataframe with all the good (value < th) file entries.
        # NOTE: after the last modifications, broken data are not present anymore.
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
    
        # build the spectrum (and spectral intensity)
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
        print('Number of photons emitted per primary (single-photon approx):', Ne, '\n')
        
        return Esteps, spectrum, Nevents, Nspectra, Nbroken, Ne
    else:
        with open(filename, 'r') as f:
            header = f.readline()
            print("reading BK file: %s.dat" % filename)
            print("header:", header)       
        data_read = np.loadtxt(filename, dtype='float', \
                               delimiter=' ', skiprows=1, unpack=True)
        Esteps = data_read[0]
        spectrum = data_read[1]
        print("file read!")       
        return Esteps, spectrum, None, None, None, None


####################### process_photons_optimized (different versions) ###########################
def process_photons_optimized(rf):
    """
    Optimized function to process photon data from ROOT file.
    
    Returns:
    - Eph: Photon energies in GeV
    - thetaX_ph: Photon angles in mrad
    - thetaY_ph: Photon angles in mrad
    - Nph: Number of photons
    """
    
    tstart = time.time()
    
    # Load only necessary columns and apply filters during loading
    branches_red = ["eventID", "volume", "angle_x", "angle_y", "Ekin", "parentID", "particle"]
    
    # Use array filtering for maximum performance
    arrays = rf['scoring_ntuple'].arrays(branches_red, library='np')
    
    # Create boolean masks using numpy (much faster than pandas filtering)
    volume_mask = arrays['volume'] == b"Detector"
    parent_mask = arrays['parentID'] > 0
    particle_mask = arrays['particle'] == b"gamma"
    
    # Combine all masks
    photon_mask = volume_mask & parent_mask & particle_mask
    
    # Apply mask and convert in one step
    Eph = arrays['Ekin'][photon_mask] * 0.001         # MeV -> GeV
    thetaX_ph = arrays['angle_x'][photon_mask] * 1e3  # rad -> mrad
    thetaY_ph = arrays['angle_y'][photon_mask] * 1e3  # rad -> mrad
    Nph = len(Eph)
    
    telapsed = time.time() - tstart
    print(f"Number of emitted photons: {Nph}")
    print(f"Elapsed time: {telapsed:.2f} s\n")
    
    return Eph, thetaX_ph, thetaY_ph, Nph


# Even faster version using uproot's built-in filtering
def process_photons_uproot_filter(rf):
    """
    Ultra-fast version using uproot's expression filtering.
    """
    
    tstart = time.time()
    
    # Use uproot's expression filtering to load only photon data
    expression = "(volume == 'Detector') & (parentID > 0) & (particle == 'gamma')"
    
    # Load only the filtered data directly
    df_ph = rf['scoring_ntuple'].arrays(
        ["Ekin", "angle_x", "angle_y"],
        cut=expression,
        library='pd'
    )
    
    # Convert to arrays and apply unit conversions
    Eph = df_ph['Ekin'].values * 0.001         # MeV -> GeV
    thetaX_ph = df_ph['angle_x'].values * 1e3  # rad -> mrad
    thetaY_ph = df_ph['angle_y'].values * 1e3  # rad -> mrad
    Nph = len(Eph)
    
    telapsed = time.time() - tstart
    print(f"Number of emitted photons: {Nph}")
    print(f"Elapsed time: {telapsed:.2f} s\n")
    
    return Eph, thetaX_ph, thetaY_ph, Nph


# Memory-efficient version for very large datasets
def process_photons_memory_efficient(rf, chunk_size=100000):
    """
    Memory-efficient version that processes data in chunks.
    """
    
    tstart = time.time()
    
    total_entries = int(rf['scoring_ntuple'].num_entries)
    Eph_list, thetaX_list, thetaY_list = [], [], []
    
    # Process data in chunks to reduce memory usage
    for start_idx in range(0, total_entries, chunk_size):
        end_idx = min(start_idx + chunk_size, total_entries)
        
        arrays = rf['scoring_ntuple'].arrays(
            ["volume", "angle_x", "angle_y", "Ekin", "parentID", "particle"],
            entry_start=start_idx,
            entry_stop=end_idx,
            library='np'
        )
        
        # Apply filters
        mask = (arrays['volume'] == b"Detector") & \
                (arrays['parentID'] > 0) & (arrays['particle'] == b"gamma")
        
        # Store results
        if np.any(mask):
            Eph_list.append(arrays['Ekin'][mask] * 0.001)
            thetaX_list.append(arrays['angle_x'][mask] * 1e3)
            thetaY_list.append(arrays['angle_y'][mask] * 1e3)
    
    # Combine all chunks
    Eph = np.concatenate(Eph_list) if Eph_list else np.array([])
    thetaX_ph = np.concatenate(thetaX_list) if thetaX_list else np.array([])
    thetaY_ph = np.concatenate(thetaY_list) if thetaY_list else np.array([])
    Nph = len(Eph)
    
    telapsed = time.time() - tstart
    print(f"Number of emitted photons: {Nph}")
    print(f"Elapsed time: {telapsed:.2f} s\n")
    
    return Eph, thetaX_ph, thetaY_ph, Nph


# One-liner version for maximum simplicity
def process_photons_quick(rf):
    """Quick and simple version for standard use cases."""
    
    tstart = time.time()
    
    # Load and filter in one step, then convert
    df_ph = rf['scoring_ntuple'].arrays(
        ["Ekin", "angle_x", "angle_y", "volume", "parentID", "particle"],
        library='pd'
    ).query('volume == "Detector" and parentID > 0 and particle == "gamma"')
    
    Eph = df_ph['Ekin'].values * 0.001
    thetaX_ph = df_ph['angle_x'].values * 1e3
    thetaY_ph = df_ph['angle_y'].values * 1e3
    Nph = len(Eph)
    
    print(f"Number of emitted photons: {Nph}")
    print(f"Elapsed time: {time.time() - tstart:.2f} s\n")
    
    return Eph, thetaX_ph, thetaY_ph, Nph
##################################################################################################

######################### load_and_process_data (different versions) #############################
def load_and_process_data(rf, Nmax=100000):
    """
    Optimized function to load and process simulation data.
    
    Parameters:
    - rf: Root file object
    - Nmax: Maximum number of events to process
    
    Returns:
    - df_in_primary_sel: Filtered input data
    - df_out_primary_sel: Filtered output data
    """
    
    # Disable warnings
    import warnings
    warnings.simplefilter("ignore")
    
    # Load only necessary columns to reduce memory
    columns_needed = ['eventID', 'angle_x', 'angle_y', 'Ekin', 'volume', 'parentID']
    df = rf['scoring_ntuple'].arrays(columns_needed, library='pd')
    
    # Create boolean masks for filtering (much faster than separate operations)
    crystal_mask = (df['volume'] == "Crystal") & (df['parentID'] == 0)
    detector_mask = (df['volume'] == "Detector") & (df['parentID'] == 0)
    
    # Apply masks and limit to Nmax in one step
    df_in_all_primary = df[crystal_mask].head(Nmax)
    df_out_all_primary = df[detector_mask].head(Nmax)
    
    # Use the actual minimum length
    actual_Nmax = min(len(df_in_all_primary), len(df_out_all_primary), int(Nmax))
    
    # Select final datasets with proper slicing
    df_in_primary_sel = df_in_all_primary[["eventID", "angle_x", "angle_y", "Ekin"]].iloc[:actual_Nmax]
    df_out_primary_sel = df_out_all_primary[["eventID", "angle_x", "angle_y", "Ekin"]].iloc[:actual_Nmax]
    
    print(f"Number of considered events (oriented case): {len(df_in_primary_sel)}")
    
    return df_in_primary_sel, df_out_primary_sel


# Even faster version using numpy arrays directly
def load_and_process_data_fast(rf, Nmax=100000):
    """
    Ultra-fast version using numpy arrays directly.
    """

    # Disable warnings
    import warnings
    warnings.simplefilter("ignore")
    
    # Load data as numpy arrays for maximum performance
    arrays = rf['scoring_ntuple'].arrays(['eventID', 'angle_x', 'angle_y', \
                                          'Ekin', 'volume', 'parentID'])
    
    # Convert to numpy operations
    volume = arrays['volume']
    parentID = arrays['parentID']
    eventID = arrays['eventID']
    angle_x = arrays['angle_x']
    angle_y = arrays['angle_y']
    Ekin = arrays['Ekin']
    
    # Create masks
    crystal_mask = (volume == b"Crystal") & (parentID == 0)
    detector_mask = (volume == b"Detector") & (parentID == 0)
    
    # Apply masks and limit to Nmax
    crystal_indices = np.where(crystal_mask)[0][:Nmax]
    detector_indices = np.where(detector_mask)[0][:Nmax]
    
    # Create final DataFrames
    df_in_primary_sel = pd.DataFrame({
        'eventID': eventID[crystal_indices],
        'angle_x': angle_x[crystal_indices],
        'angle_y': angle_y[crystal_indices],
        'Ekin': Ekin[crystal_indices]
    })
    
    df_out_primary_sel = pd.DataFrame({
        'eventID': eventID[detector_indices],
        'angle_x': angle_x[detector_indices],
        'angle_y': angle_y[detector_indices],
        'Ekin': Ekin[detector_indices]
    })
    
    print(f"Number of considered events (oriented case): {len(df_in_primary_sel)}")
    
    return df_in_primary_sel, df_out_primary_sel


# Memory-optimized version for very large files
def load_and_process_data_memory_optimized(rf, Nmax=100000, chunk_size=50000):
    """
    Memory-optimized version that processes data in chunks.
    Useful for very large files that don't fit in memory.
    """
    warnings.simplefilter("ignore")
    
    # Get the total number of entries to pre-allocate arrays
    n_entries = rf['scoring_ntuple'].num_entries
    
    # Pre-allocate lists (more memory efficient than growing arrays)
    in_events = []
    out_events = []
    
    # Process in chunks to reduce memory usage
    for start in range(0, n_entries, chunk_size):
        end = min(start + chunk_size, n_entries)
        
        df_chunk = rf['scoring_ntuple'].arrays(
            ['eventID', 'angle_x', 'angle_y', 'Ekin', 'volume', 'parentID'],
            entry_start=start,
            entry_stop=end,
            library='pd'
        )
        
        # Filter and collect events
        crystal_events = df_chunk[(df_chunk['volume'] == "Crystal") & (df_chunk['parentID'] == 0)]
        detector_events = df_chunk[(df_chunk['volume'] == "Detector") & (df_chunk['parentID'] == 0)]
        
        in_events.append(crystal_events[["eventID", "angle_x", "angle_y", "Ekin"]])
        out_events.append(detector_events[["eventID", "angle_x", "angle_y", "Ekin"]])
        
        # Stop if we have enough events
        if len(in_events) > 0 and sum(len(chunk) for chunk in in_events) >= Nmax:
            break
    
    # Concatenate and limit to Nmax
    df_in_primary_sel = pd.concat(in_events, ignore_index=True).head(Nmax)
    df_out_primary_sel = pd.concat(out_events, ignore_index=True).head(Nmax)
    
    print(f"Number of considered events (oriented case): {len(df_in_primary_sel)}")
    
    return df_in_primary_sel, df_out_primary_sel
##################################################################################################

######################### get deflection angles (different versions) #############################
def get_deflection_angles_old(df_out, df_in, pot_good_events, 
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


def get_deflection_angles(df_merged, 
                          res_thetaX_in=0, res_thetaY_in=0, \
                          res_DthetaX=0, res_DthetaY=0, \
                          ang_cut=np.pi*1e6, cut_center=[0., 0.]):
    """
    Calculate the deflection for events whose transverse
    deviation angle is within a given angular cut.
    Experimental resolution can be taken into account.
    Input angles must be given in [urad]. However, the angles
    in the passed merged dataframe are expressed in [rad].
    In this case, there is no neee dof pot_good_events list.
    """
    pot_good_events = df_merged.eventID.unique()
    DthetaX = [] #urad
    DthetaY = [] #urad
    thetaX_in = [] #rad
    thetaY_in = [] #rad
    good_events = []
    df_grouped = df_merged.groupby('eventID')
    for i, df_i in tqdm(df_grouped):
        thX_in = np.random.normal(df_i['angle_x_in'].values[0]*1e6, res_thetaX_in, 1)[0] #urad
        thY_in = np.random.normal(df_i['angle_y_in'].values[0]*1e6, res_thetaY_in, 1)[0] #urad
        if np.sqrt((thX_in-cut_center[0])**2 + (thY_in-cut_center[1])**2) < ang_cut: #urad
            thetaX_in.append(thX_in)
            thetaY_in.append(thY_in)
            DthetaX_mean = (df_i['angle_x_out'].values[0] - df_i['angle_x_in'].values[0])*1e6 #urad
            DthetaY_mean = (df_i['angle_y_out'].values[0] - df_i['angle_y_in'].values[0])*1e6 #urad  
            DthetaX.append(np.random.normal(DthetaX_mean, res_DthetaX, 1)[0])
            DthetaY.append(np.random.normal(DthetaY_mean, res_DthetaY, 1)[0])
            good_events.append(i)
        #del df_i
    return DthetaX, DthetaY, thetaX_in, thetaY_in, good_events


def get_deflection_angles_vectorized(df_merged, 
                                     res_thetaX_in=0, 
                                     res_thetaY_in=0,
                                     res_DthetaX=0, 
                                     res_DthetaY=0,
                                     ang_cut=np.pi*1e6, 
                                     cut_center=(0., 0.)):
    """
    Calculate deflection angles for events within angular cut.
    Vectorized implementation for much better performance.
    
    Parameters:
    - df_merged: DataFrame with eventID, angle_x_in, angle_y_in, angle_x_out, angle_y_out
    - res_*: Resolution uncertainties (standard deviation) in urad
    - ang_cut: Angular cut radius in urad
    - cut_center: Center of angular cut in urad
    
    Returns:
    - DthetaX, DthetaY: Deflection angles with resolution smearing
    - thetaX_in, thetaY_in: Input angles with resolution smearing  
    - good_events: Event IDs that pass the cut
    """
    # Convert to urad once (vectorized)
    theta_x_in_urad = df_merged['angle_x_in'].values * 1e6
    theta_y_in_urad = df_merged['angle_y_in'].values * 1e6
    
    # Calculate deflection in urad
    DthetaX_mean = (df_merged['angle_x_out'].values - df_merged['angle_x_in'].values) * 1e6
    DthetaY_mean = (df_merged['angle_y_out'].values - df_merged['angle_y_in'].values) * 1e6
    
    # Apply resolution smearing
    if res_thetaX_in > 0 or res_thetaY_in > 0:
        theta_x_in_smeared = np.random.normal(theta_x_in_urad, res_thetaX_in)
        theta_y_in_smeared = np.random.normal(theta_y_in_urad, res_thetaY_in)
    else:
        theta_x_in_smeared = theta_x_in_urad
        theta_y_in_smeared = theta_y_in_urad
    
    # Calculate radial distance from cut center
    radial_dist = np.sqrt((theta_x_in_smeared - cut_center[0])**2 + 
                          (theta_y_in_smeared - cut_center[1])**2)
    
    # Apply angular cut
    mask = radial_dist < ang_cut
    good_events = df_merged['eventID'].values[mask]
    
    # Apply resolution to deflection angles
    if res_DthetaX > 0 or res_DthetaY > 0:
        DthetaX = np.random.normal(DthetaX_mean[mask], res_DthetaX)
        DthetaY = np.random.normal(DthetaY_mean[mask], res_DthetaY)
    else:
        DthetaX = DthetaX_mean[mask]
        DthetaY = DthetaY_mean[mask]
    
    thetaX_in = theta_x_in_smeared[mask]
    thetaY_in = theta_y_in_smeared[mask]
    
    return DthetaX, DthetaY, thetaX_in, thetaY_in, good_events


# Alternative version for even better performance with pre-grouped data
def get_deflection_angles_optimized(df_merged,
                                    res_thetaX_in=0,
                                    res_thetaY_in=0, 
                                    res_DthetaX=0,
                                    res_DthetaY=0,
                                    ang_cut=np.pi*1e6,
                                    cut_center=(0., 0.)):
    """
    Optimized version that assumes df_merged has one row per eventID.
    If there are duplicates, use df_merged = df_merged.drop_duplicates('eventID') first.
    """
    # Convert to urad
    theta_x_in_urad = df_merged['angle_x_in'] * 1e6
    theta_y_in_urad = df_merged['angle_y_in'] * 1e6
    
    # Calculate deflections
    DthetaX_mean = (df_merged['angle_x_out'] - df_merged['angle_x_in']) * 1e6
    DthetaY_mean = (df_merged['angle_y_out'] - df_merged['angle_y_in']) * 1e6
    
    # Apply resolution smearing to input angles
    if res_thetaX_in > 0:
        theta_x_in_urad = np.random.normal(theta_x_in_urad, res_thetaX_in)
    if res_thetaY_in > 0:
        theta_y_in_urad = np.random.normal(theta_y_in_urad, res_thetaY_in)
    
    # Apply angular cut
    radial_dist = np.hypot(theta_x_in_urad - cut_center[0], 
                           theta_y_in_urad - cut_center[1])
    mask = radial_dist < ang_cut
    
    # Apply resolution to deflection angles
    if res_DthetaX > 0:
        DthetaX = np.random.normal(DthetaX_mean[mask], res_DthetaX)
    else:
        DthetaX = DthetaX_mean[mask]
        
    if res_DthetaY > 0:
        DthetaY = np.random.normal(DthetaY_mean[mask], res_DthetaY)
    else:
        DthetaY = DthetaY_mean[mask]
    
    return (DthetaX, DthetaY, 
            theta_x_in_urad[mask], theta_y_in_urad[mask], 
            df_merged['eventID'].values[mask])


# If you need to handle multiple entries per eventID, use this first:
def prepare_dataframe(df_merged):
    """Ensure one row per eventID by taking first occurrence"""
    return df_merged.groupby('eventID').first().reset_index()
##################################################################################################

######################### filter photon emission using deflection values #########################
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
    I updated the function to consider also 2 or 3 circular selctions, as in the Bent Ge<110> article. 
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
    Erad_sel = np.array([sum(Eph_dict_sel[event]) for event in Eph_dict_sel.keys()]) 
    #Erad_Sel is the total energy lost by radiation per (selected) event [GeV]
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
############# Set of functions useful for TestBeamOC(2) ###############################################
#######################################################################################################
def applyCutsToDF(df, th_APC1=1e15, th_APC2=0, crystal_width_cm=1e15, crystal_height_cm=1e15, \
                  angular_cut_type='theta', theta_cut=3.1415):
    """
    It applies cuts to a df that contains at least x, y, thetax, thetay, APC1, APC2.
    The following variables are a subset of the columns of outData ntuple of TestBeamOC(2).
    The angular cut can be applied on thetaX, thetaY or on theta, which is valid in the axial case.
    """
    import numpy as np
    
    if angular_cut_type == "thetaX":
        ang_cut = np.abs(df.thetaX) < theta_cut
    elif angular_cut_type == "thetaY":
        ang_cut = np.abs(df.thetaY) < theta_cut
    else:
        ang_cut = df.thetaX**2 + df.thetaY**2 < theta_cut**2
    print("applied %.2e rad angular cut on %s" % (theta_cut, angular_cut_type))
           
    if 'edep_APC1' not in df.columns:
        df_sel = df[(np.abs(df.x_cry) < crystal_width_cm*0.5) & \
                    (np.abs(df.y_cry) < crystal_height_cm*0.5) & (ang_cut)]
    else:
        if len(df[df.edep_APC1 == 0]) == len(df.edep_APC1):
            df_sel = df[(np.abs(df.x_cry) < crystal_width_cm*0.5) & \
                        (np.abs(df.y_cry) < crystal_height_cm*0.5) & (ang_cut)]
        else:
            df_sel = df[(df.edep_APC1 < th_APC1) & (df.edep_APC2 > th_APC2) & \
                        (np.abs(df.x_cry) < crystal_width_cm*0.5) & \
                        (np.abs(df.y_cry) < crystal_height_cm*0.5) & (ang_cut)]
    return df_sel


def FindPartPerEvent(df):
    """
    It counts the number of particles per event present
    in an ntuple that contains 'eventID' column, 
    such as 'scoringScreen' of TestBeamOC.
    WARNING: there may be particles that impinge multiple times 
    on the screen (due to scattering), thus it is an approx.
    """
    import pandas as pd
    import numpy as np
    from tqdm.notebook import tqdm
    part_events = []
    df_dict = {}
    df_grouped = df.groupby('eventID')
    for i, df_i in tqdm(df_grouped):
        part_events.append(len(df_i))
        df_dict[i] = df_i    
    return np.array(part_events), df_dict


def CountParticles(df):
    """
    It counts the number of particles of each type present
    in an ntuple that contains 'particle' column, 
    such as 'scoringScreen' of TestBeamOC.
    WARNING: there may be particles that impinge multiple times 
    on the screen (due to scattering), thus it is an approx.
    """
    import pandas as pd
    import numpy as np
    from tqdm.notebook import tqdm
    part_counts = {}
    df_grouped = df.groupby('particle')
    for i, df_i in tqdm(df_grouped):
        part_counts[i] = len(df_i)
    return part_counts
#######################################################################################################

#######################################################################################################
#### Functions implemented by R. Negrello to fit theta_out curve considering various contributions ####
#######################################################################################################
def f_dech(theta_X, theta_dech, A_dech, path_file_json="parameters.json"):
    if os.path.exists(path_file_json):
        with open(path_file_json, 'r') as file:
            parameters = json.load(file)
        sigma_VR = parameters["sigma"]  
        theta_VR = parameters["mu"]
        theta_Ch = parameters["mu"]
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
    else:
        return 0
    
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
#######################################################################################################
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
    To avoid future ambiguity, I changed the second column name from 'state' to
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
#######################################################################################################
def merge_root_ntuples(data_path, ntuple):
    """
    This function is quite general, it merges ntuples in python.
    There are no constranints on the filenames nor on the eventID.
    It returns a dataframe.
    """
    import uproot
    dataframes = []
    for filename in os.listdir(data_path):
        print('opening', filename)
        file_path = os.path.join(data_path, filename)
        if os.path.isfile(file_path):
            rf = uproot.open(file_path)
            rf_content = [item.split(';')[0] for item in rf.keys()]
            if ntuple in rf_content:
                df = rf[ntuple_name].arrays(library='pd')
            dataframes.append(df)
    df_merged = pd.concat(dataframes, ignore_index=True)
    print('ntuples merged into a single dataframe!')
    return df_merged


def merge_root_ntuples_based_on_eventID(data_path, ntuple_name, beVerbose=False):
    """
    Function to merge an ntuple contained in many root files.
    During the merging, it automatically corrects the eventID, 
    based on the number of events progressively read.
    It returns a dataframe with the merged ntuple.
    WARNING: if you have a file with more than one ntuple, 
    by reading them separately, you'll lose the correlation.
    OLD (full gpaterno) VERSION.
    """
    import uproot
    dataframes_list = []
    add2eventIDc = 0
    add2eventIDi = 0
    events_read = 0
    nfiles = 0
    add1 = 0
    fIDs = [int(filename.split('.root')[0].split('_')[-1]) for filename in os.listdir(data_path)]
    if beVerbose:
        print("ID of files to merge:", fIDs, "\n")
    for filename in os.listdir(data_path):
        if beVerbose:
            print('opening', filename, '...')    
        # open a root file
        file = os.path.join(data_path, filename)
        rf = uproot.open(file)
        rf_content = [item.split(';')[0] for item in rf.keys()]      
        # retrieve the ntuple
        df_root = rf[ntuple_name].arrays(library='pd')
        # list of events        
        event_list = list(df_root.eventID.unique())
        #print("event_list:", event_list)
        n_events = len(event_list)
        # correct the event number and create a dataframe
        add2eventIDi = max([n_events, max(event_list)])
        for event in event_list:
            df_root.loc[df_root.eventID == event, 'eventID'] += add2eventIDc + add1 if min(event_list)==0 else 0
        add1 = 1 if (n_events < max(event_list)) else 0
        #print("add1:", add1)
        dataframes_list.append(df_root)          
        # increase te number of files read
        add2eventIDc += add2eventIDi
        events_read += n_events
        nfiles += 1
        if beVerbose:
            print("n_events:", n_events)
            print("add2eventIDc:", add2eventIDc)
    # merge the dataframes
    df_root_merged = pd.concat(dataframes_list, ignore_index=True)
    # return the merged dataframe
    print("\nevents_read:", events_read)
    print('%d files merged!\n' % (nfiles))
    return df_root_merged 


def merge_FastSimChannelingRad_files_old(data_path, correct_particle=False, primary='e-', beVerbose=False, save_result=False):
    """
    Function to merge both the root and the txt files of a set of Geant4
    simulations of particle interaction in Oriented Crystals, obatined with
    FastSimChannelingRad application.
    During the merging, it automatically corrects the eventID.
    The name of the corresponding root and txt files MUST be the same.
    It returns two dataframes (one for each ntuple) with the merged data, 
    but it can also create a merged root file.
    NOTE1: due to some issues in the creation of the output tree with
    branches of strings, I converted strings in integers
    (look at volume_dict and part_dict to know the coding).
    OLD (full gpaterno) VERSION.
    """
    import uproot
    dataframes_defl_list = []
    dataframes_rad_list = []
    dataframes_txt_list = []
    dataframes_edep_calo_list = []
    dataframes_edep_vol_list = []
    add2eventIDc = 0
    add2eventIDi = 0
    events_read = 0
    nfiles = 0
    add1 = 0
    file_list = os.listdir(data_path)
    fIDs = [int(filename.split('.root')[0].split('_')[-1]) for filename in file_list]
    if beVerbose:
        print("ID of files to merge:", fIDs, "\n")
    for filename in file_list:
        if beVerbose:
            print('opening', filename, '...')    
        # open a root file
        file = os.path.join(data_path, filename)
        rf = uproot.open(file)
        rf_content = [item.split(';')[0] for item in rf.keys()]      
        # deflection ntuple
        df_defl = rf['scoring_ntuple'].arrays(library='pd')          
        # radiation ntuple
        df_rad = rf['photon_spectrum'].arrays(library='pd')
        #df_rad = df_rad.rename(columns={"E": "Ekin"})
        # edep_vol calo
        if 'edep_calo' in rf_content:
            df_edep_calo = rf['edep_calo'].arrays(library='pd')
        else:
             df_edep_calo = pd.DataFrame()
        # edep_vol ntuple
        if 'edep_vol' in rf_content:
            df_edep_vol = rf['edep_vol'].arrays(library='pd')    
        else:
             df_edep_vol = pd.DataFrame()
        # list of events
        event_list = list(df_defl.eventID.unique())
        #print("event_list:", event_list)
        n_events = len(event_list)
        # correct the event number and create a dataframe
        add2eventIDi = max([n_events, max(event_list)])
        for event in event_list:
            df_defl.loc[df_defl.eventID == event, 'eventID'] += add2eventIDc + add1 if min(event_list)==0 else 0
            if not df_rad.empty:
                df_rad.loc[df_rad.eventID == event, 'eventID'] += add2eventIDc + add1 if min(event_list)==0 else 0
            if not df_edep_calo.empty:
                df_edep_calo.loc[df_edep_calo.eventID == event, 'eventID'] += add2eventIDc + add1 if min(event_list)==0 else 0
            if not df_edep_vol.empty:
                df_edep_vol.loc[df_edep_vol.eventID == event, 'eventID'] += add2eventIDc + add1 if min(event_list)==0 else 0
        add1 = 1 if (n_events < max(event_list)) else 0
        #print("add1:", add1)
        dataframes_defl_list.append(df_defl)
        if not df_rad.empty:
            dataframes_rad_list.append(df_rad)
        if not df_edep_calo.empty:
            dataframes_edep_calo_list.append(df_edep_calo)
        if not df_edep_vol.empty:
            dataframes_edep_vol_list.append(df_edep_vol) 
        # open txt file
        volume = []
        eventID = []
        trackID = [] 
        x = []
        tx = []
        z = []
        xx = []
        yy = []
        file_tag = os.path.join(data_path.replace('rad', 'tag'), filename.replace('root', 'txt'))
        try:
            with open(file_tag, 'r') as txt_file:
                if beVerbose:
                    #print('opening', file_tag, '...')
                    print('opening', filename.replace('root', 'txt'), '...')
                for line in txt_file:
                    values = line.strip().split()
                    eventID.append(int(values[0]) + events_read)
                    trackID.append(int(values[1]))
                    x.append(float(values[2]))
                    tx.append(float(values[3]))
                    z.append(float(values[4]))
                    xx.append(float(values[5]))
                    if len(values) > 6:
                        yy.append(float(values[6]))
                    else:
                        yy.append(0)
            df = pd.DataFrame({
                "volume": "Crystal",
                "eventID": eventID,
                "trackID": trackID,
                "x": x,
                "tx": tx,
                "z": z,
                "xx": xx,
                "yy": yy
            })
            dataframes_txt_list.append(df)  
        except:
            if beVerbose:
                print('file %s not found!' % filename.replace('root', 'txt'))
            else:
                continue      
        # increase te number of files read
        add2eventIDc += add2eventIDi
        events_read += n_events
        nfiles += 1
        if beVerbose:
            print("n_events:", n_events)
            print("add2eventIDc:", add2eventIDc)
    # merge the dataframes
    df_defl_merged = pd.concat(dataframes_defl_list, ignore_index=True)
    if correct_particle:
        df_defl_merged.particle = [primary]*len(df_defl_merged) #correction #<-- it may not be required!
    df_rad_merged = pd.concat(dataframes_rad_list, ignore_index=True)
    if not len(dataframes_edep_calo_list) == 0:
        df_edep_calo_merged = pd.concat(dataframes_edep_calo_list, ignore_index=True)
    else:
        df_edep_calo_merged = pd.DataFrame({'eventID' : [], 'edep' : []})
    if not len(dataframes_edep_vol_list) == 0:
        df_edep_vol_merged = pd.concat(dataframes_edep_vol_list, ignore_index=True)
    else:
        df_edep_vol_merged = pd.DataFrame({'eventID' : [], 'volID' : [] , 'edep' : []})
    if not len(dataframes_txt_list) == 0:
        df_txt_merged = pd.concat(dataframes_txt_list, ignore_index=True)
    else:
        df_txt_merged = pd.DataFrame({'eventID' : []})
    # save merged dafatrames to proper files
    if save_result:
        # export merged defl and rad dataframes to root files
        last_file = file_list[-1].replace(".root", "")
        if '_' in last_file:
            outputfile = last_file[:last_file.find("_", -1)] + 'merged' + str(nfiles ) + 'files'
        else:
            outputfile = 'merged' + str(nfiles ) + 'files'
        rf_merged = uproot.recreate(data_path + '../' + outputfile + '.root')
        tree_defl = rf_merged.mktree("scoring_ntuple", {
                                                 "eventID": np.int32, "volume": np.int32, \
                                                 "x": np.float64, "y": np.float64, \
                                                 "angle_x": np.float64, "angle_y": np.float64, \
                                                 "Ekin": np.float64, "particle": np.int32, \
                                                 "particleID": np.int32, "parentID": np.int32
                                                })
        # conversion of string columns into integer columns (without this an error occurs)
        df_defl_merged["volume"] = df_defl_merged["volume"].astype(str)
        volume_dict = {item: i for i,item in enumerate(df_defl_merged.volume.unique())}
        df_defl_merged["volume"] = df_defl_merged["volume"].replace(volume_dict)
        df_defl_merged["particle"] = df_defl_merged["particle"].astype(str) 
        part_dict = {item: i for i,item in enumerate(df_defl_merged.particle.unique())}
        df_defl_merged["particle"] = df_defl_merged["particle"].replace(part_dict)
        #tree_defl
        tree_defl.extend({
                          "eventID": df_defl_merged.eventID.values, "volume": df_defl_merged.volume.values, \
                          "x": df_defl_merged.x.values, "y": df_defl_merged.y.values, \
                          "angle_x": df_defl_merged.angle_x.values, "angle_y": df_defl_merged.angle_y.values, 
                          "Ekin": df_defl_merged.Ekin.values, "particle": df_defl_merged.particle.values, \
                          "particleID": df_defl_merged.particleID.values, "parentID": df_defl_merged.parentID.values
                         })        
        #tree_rad
        tree_rad = rf_merged.mktree("photon_spectrum", {
                                                 "E": np.float64, \
                                                 "angle_x": np.float64, "angle_y": np.float64, \
                                                 "eventID": np.int32,
                                                })    
        tree_rad.extend({
                         "E": df_rad_merged.E.values, \
                         "angle_x": df_rad_merged.angle_x.values, \
                         "angle_y": df_rad_merged.angle_y.values, \
                         "eventID": df_rad_merged.eventID.values
                        })
        #tree_edep_calo
        tree_edep_calo = rf_merged.mktree("edep_calo", {
                                                 "eventID": np.int32, \
                                                 "edep": np.float64
                                                }) 
        tree_edep_calo.extend({
                         "eventID": df_edep_calo_merged.eventID.values, \
                         "edep": df_edep_calo_merged.edep
                        })
        #tree_edep_vol
        tree_edep_vol = rf_merged.mktree("edep_vol", {
                                                 "eventID": np.int32, \
                                                 "volID": np.int32, \
                                                 "edep": np.float64
                                                }) 
        tree_edep_vol.extend({
                         "eventID": df_edep_vol_merged.eventID.values, \
                         "volID": df_edep_vol_merged.volID, \
                         "edep": df_edep_vol_merged.edep
                        })
        # export merged txt dataframe to text file
        if not df_txt_merged.empty:
            with open(data_path + '../' + outputfile + '.txt', 'a') as f:
                df_string = df_txt_merged.to_string(header=False, index=False)
                f.write(df_string)
    # return the merged dataframes
    print("\n")
    if save_result:
        print("volume_dict merged:", volume_dict)
        print("part_dict merged:", part_dict)
    else:
        volume_dict = {}
        part_dict = {}
    print("events_read:", events_read)
    print('%d files merged!\n' % (nfiles))
    return df_defl_merged, df_rad_merged, df_txt_merged, df_edep_calo_merged, df_edep_vol_merged, volume_dict, part_dict


def merge_FastSimChannelingRad_files(data_path, correct_particle=False, primary='e-', beVerbose=False, save_result=False):
    """
    Optimized (DeepSeek) function to merge root files of Geant4 simulations.
    """
    
    import uproot
    import numpy as np
    import pandas as pd
    import os
    
    # Pre-allocate lists for dataframes (for faster appending)
    dataframes_defl_list = []
    dataframes_rad_list = []
    dataframes_txt_list = []
    dataframes_edep_calo_list = []
    dataframes_edep_vol_list = []
    
    events_read = 0
    nfiles = 0
    
    # Get file list and sort to ensure consistent processing
    file_list = sorted([f for f in os.listdir(data_path) if f.endswith('.root')])
    
    if beVerbose:
        print(f"Number of files to merge: {len(file_list)}")
    
    # Precompute file paths and IDs (vectorized)
    root_files = [os.path.join(data_path, f) for f in file_list]
    txt_files = [os.path.join(data_path.replace('rad', 'tag'), f.replace('.root', '.txt')) for f in file_list]
    
    # Process files
    for root_file, txt_file, filename in zip(root_files, txt_files, file_list):
        if beVerbose:
            print(f'Processing {filename}...')
        
        # Open ROOT file
        try:
            rf = uproot.open(root_file)
        except:
            if beVerbose:
                print(f"Failed to open {root_file}")
            continue
        
        # Load data in bulk - single operation per file
        # Process deflection ntuple
        df_defl = rf['scoring_ntuple'].arrays(library='pd')
        n_events = df_defl['eventID'].max() + 1  # Assuming eventID starts at 0
        
        # Update eventIDs in bulk (vectorized)
        df_defl['eventID'] += events_read
        
        # Process radiation ntuple if exists
        df_rad = pd.DataFrame()
        if 'photon_spectrum' in rf:
            df_rad = rf['photon_spectrum'].arrays(library='pd')
            #df_rad = df_rad.rename(columns={"E": "Ekin"})
            if not df_rad.empty:
                df_rad['eventID'] += events_read
        
        # Process edep_calo if exists
        df_edep_calo = pd.DataFrame()
        if 'edep_calo' in rf:
            df_edep_calo = rf['edep_calo'].arrays(library='pd')
            if not df_edep_calo.empty:
                df_edep_calo['eventID'] += events_read
        
        # Process edep_vol if exists
        df_edep_vol = pd.DataFrame()
        if 'edep_vol' in rf:
            df_edep_vol = rf['edep_vol'].arrays(library='pd')
            if not df_edep_vol.empty:
                df_edep_vol['eventID'] += events_read
        
        # Append dataframes (store references, not copies)
        dataframes_defl_list.append(df_defl)
        if not df_rad.empty:
            dataframes_rad_list.append(df_rad)
        if not df_edep_calo.empty:
            dataframes_edep_calo_list.append(df_edep_calo)
        if not df_edep_vol.empty:
            dataframes_edep_vol_list.append(df_edep_vol)
        
        # Process text file if exists
        if os.path.exists(txt_file):
            if beVerbose:
                print(f'Processing text file {os.path.basename(txt_file)}...')
            
            # Read text file using numpy for speed
            try:
                data = np.loadtxt(txt_file)
                if len(data.shape) == 1:
                    data = data.reshape(1, -1)
                
                # Create dataframe from numpy array
                n_cols = data.shape[1]
                if n_cols >= 7:
                    df_txt = pd.DataFrame({
                        "volume": "Crystal",
                        "eventID": data[:, 0].astype(int) + events_read,
                        "trackID": data[:, 1].astype(int),
                        "x": data[:, 2].astype(float),
                        "tx": data[:, 3].astype(float),
                        "z": data[:, 4].astype(float),
                        "xx": data[:, 5].astype(float),
                        "yy": data[:, 6].astype(float) if n_cols > 6 else np.zeros(len(data))
                    })
                elif n_cols == 6:
                    df_txt = pd.DataFrame({
                        "volume": "Crystal",
                        "eventID": data[:, 0].astype(int) + events_read,
                        "trackID": data[:, 1].astype(int),
                        "x": data[:, 2].astype(float),
                        "tx": data[:, 3].astype(float),
                        "z": data[:, 4].astype(float),
                        "xx": data[:, 5].astype(float),
                        "yy": np.zeros(len(data))
                    })
                else:
                    if beVerbose:
                        print(f"Unexpected format in {txt_file}")
                    continue
                
                dataframes_txt_list.append(df_txt)
                
            except Exception as e:
                if beVerbose:
                    print(f'Error reading {txt_file}: {e}')
        elif beVerbose:
            print(f'Text file {os.path.basename(txt_file)} not found')
        
        # Update counters
        events_read += n_events
        nfiles += 1
        
        if beVerbose:
            print(f"Events in this file: {n_events}")
            print(f"Total events read: {events_read}")
    
    # Merge dataframes (use concat only once)
    if dataframes_defl_list:
        df_defl_merged = pd.concat(dataframes_defl_list, ignore_index=True, copy=False)
    else:
        df_defl_merged = pd.DataFrame()
    
    if correct_particle:
        df_defl_merged['particle'] = primary
    
    df_rad_merged = pd.concat(dataframes_rad_list, ignore_index=True, copy=False) if dataframes_rad_list else pd.DataFrame()
    df_edep_calo_merged = pd.concat(dataframes_edep_calo_list, ignore_index=True, copy=False) if dataframes_edep_calo_list else pd.DataFrame()
    df_edep_vol_merged = pd.concat(dataframes_edep_vol_list, ignore_index=True, copy=False) if dataframes_edep_vol_list else pd.DataFrame()
    df_txt_merged = pd.concat(dataframes_txt_list, ignore_index=True, copy=False) if dataframes_txt_list else pd.DataFrame()
    
    # Save results if requested
    if save_result and nfiles > 0:
        # Generate output filename
        last_file = file_list[-1].replace(".root", "")
        if '_' in last_file:
            last_underscore = last_file.rfind('_')
            outputfile = last_file[:last_underscore] + '_merged' + str(nfiles) + 'files'
        else:
            outputfile = 'merged' + str(nfiles) + 'files'
        
        output_path = os.path.join(os.path.dirname(data_path), '../', outputfile + '.root')
        
        # Create output file
        rf_merged = uproot.recreate(output_path)
        
        # Save deflection ntuple
        if not df_defl_merged.empty:
            # Store original mappings            
            volume_dict = {item: i for i,item in enumerate(df_defl_merged.volume.unique())}
            part_dict = {item: i for i,item in enumerate(df_defl_merged.particle.unique())}

            # Convert string columns to categorical codes (faster than dict replacement)
            df_defl_merged['volume'] = pd.Categorical(df_defl_merged['volume']).codes
            df_defl_merged['particle'] = pd.Categorical(df_defl_merged['particle']).codes
            
            tree_defl = rf_merged.mktree("scoring_ntuple", {
                "eventID": np.int32, "volume": np.int32,
                "x": np.float64, "y": np.float64,
                "angle_x": np.float64, "angle_y": np.float64,
                "Ekin": np.float64, "particle": np.int32,
                "particleID": np.int32, "parentID": np.int32
            })
            
            tree_defl.extend({
                "eventID": df_defl_merged.eventID.values,
                "volume": df_defl_merged.volume.values,
                "x": df_defl_merged.x.values,
                "y": df_defl_merged.y.values,
                "angle_x": df_defl_merged.angle_x.values,
                "angle_y": df_defl_merged.angle_y.values,
                "Ekin": df_defl_merged.Ekin.values,
                "particle": df_defl_merged.particle.values,
                "particleID": df_defl_merged.particleID.values,
                "parentID": df_defl_merged.parentID.values
            })
        
        # Save radiation ntuple
        if not df_rad_merged.empty:
            tree_rad = rf_merged.mktree("photon_spectrum", {
                "E": np.float64,
                "angle_x": np.float64, 
                "angle_y": np.float64,
                "eventID": np.int32,
            })
            
            tree_rad.extend({
                "E": df_rad_merged.E.values,
                "angle_x": df_rad_merged.angle_x.values,
                "angle_y": df_rad_merged.angle_y.values,
                "eventID": df_rad_merged.eventID.values
            })
        
        # Save text data
        if not df_txt_merged.empty:
            txt_output_path = os.path.join(os.path.dirname(data_path), '../', outputfile + '.txt')
            # Save as CSV with header for easier reading
            df_txt_merged.to_csv(txt_output_path, sep=' ', index=False, header=False)
            
        # Save edep_calo ntuple
        if not df_edep_calo_merged.empty:
            tree_edep_calo = rf_merged.mktree("edep_calo", {
                                                 "eventID": np.int32,
                                                 "edep": np.float64
                                                }) 
            tree_edep_calo.extend({
                             "eventID": df_edep_calo_merged.eventID,
                             "edep": df_edep_calo_merged.edep
                            })
            
        # Save edep_vol ntuple
        if not df_edep_vol_merged.empty:
            tree_edep_vol = rf_merged.mktree("edep_vol", {
                                                     "eventID": np.int32,
                                                     "volID": np.int32, 
                                                     "edep": np.float64
                                                    }) 
            tree_edep_vol.extend({
                             "eventID": df_edep_vol_merged.eventID,
                             "volID": df_edep_vol_merged.volID, 
                             "edep": df_edep_vol_merged.edep
                            })
  
    # Print summary and the merged dataframes
    print("\n")
    if save_result:
        print("volume_dict merged:", volume_dict)
        print("part_dict merged:", part_dict)
    else:
        volume_dict = {}
        part_dict = {}
    print("events_read:", events_read)
    print('%d files merged!\n' % (nfiles))
    return df_defl_merged, df_rad_merged, df_txt_merged, df_edep_calo_merged, df_edep_vol_merged, volume_dict, part_dict


def merge_TestBeamOC_files_old(data_path, beVerbose=False, save_result=False):
    """
    Function to merge the root files of a set of Geant4 simulations 
    of particle interaction in Oriented Crystals, obatined with
    TestBeamOC or TestBeamPS applications.
    During the merging, it automatically corrects the eventID.
    It returns two dataframes (one for each ntuple) with the merged data, 
    but it can also create a merged root file.
    NOTE1: due to some issues in the creation of the output tree with
    branches of strings, I converted strings in integers
    (look at part_dict to know the coding).
    OLD (full gpaterno) VERSION.
    """
    import uproot
    dataframes_out_list = []
    dataframes_scr_list = []
    add2eventIDc = 0
    add2eventIDi = 0
    events_read = 0
    nfiles = 0
    add1 = 0
    file_list = os.listdir(data_path)
    fIDs = [int(filename.split('.root')[0].split('_')[-1]) for filename in file_list]
    if beVerbose:
        print("ID of files to merge:", fIDs, "\n")
    for filename in file_list:
        if beVerbose:
            print('opening', filename, '...')    
        # open a root file
        file = os.path.join(data_path, filename)
        rf = uproot.open(file)
        rf_content = [item.split(';')[0] for item in rf.keys()]      
        # outData ntuple
        df_out = rf['outData'].arrays(library='pd')
        # scoringScreen ntuple
        df_scr = rf['scoringScreen'].arrays(library='pd')
        # list of events        
        event_list = list(df_scr.eventID.unique())
        #print("event_list:", event_list)
        n_events = len(event_list)
        # correct the event number and create a dataframe
        add2eventIDi = max([n_events, max(event_list)])
        for event in event_list:
            df_out.loc[df_out.eventID == event, 'eventID'] += add2eventIDc + add1 if min(event_list)==0 else 0
            df_scr.loc[df_scr.eventID == event, 'eventID'] += add2eventIDc + add1 if min(event_list)==0 else 0
        add1 = 1 if (n_events < max(event_list)) else 0
        #print("add1:", add1)
        dataframes_out_list.append(df_out)
        dataframes_scr_list.append(df_scr)            
        # increase te number of files read
        add2eventIDc += add2eventIDi
        events_read += n_events
        nfiles += 1
        if beVerbose:
            print("n_events:", n_events)
            print("add2eventIDc:", add2eventIDc)
    # merge the dataframes
    df_out_merged = pd.concat(dataframes_out_list, ignore_index=True)
    df_scr_merged = pd.concat(dataframes_scr_list, ignore_index=True)
    # save merged dafatrames to proper files
    if save_result:
        # export merged defl and rad dataframes to root files
        last_file = file_list[-1].replace(".root", "")
        if '_' in last_file:
            outputfile = last_file[:last_file.find("_", -1)] + 'merged' + str(nfiles ) + 'files'
        else:
            outputfile = 'merged' + str(nfiles ) + 'files'
        rf_merged = uproot.recreate(data_path + '../' + outputfile + '.root')
        tree_out = rf_merged.mktree("outData", {
                                                "eventID": np.int32, \
                                                "Tracker_NHit_X_1": np.int32, "Tracker_NHit_Y_1": np.int32, \
                                                "Tracker_NHit_X_2": np.int32, "Tracker_NHit_Y_2": np.int32, \
                                                "Tracker_X_1": np.float64, "Tracker_Y_1": np.float64, \
                                                "Tracker_X_2": np.float64, "Tracker_Y_2": np.float64, \
                                                "Ekin": np.float64, \
                                                "edep_APC1": np.float64, "edep_APC2": np.float64, \
                                                "edep_calo": np.float64, "edep_calo2": np.float64, \
                                                "edep_calo3": np.float64, "edep_screen": np.float64, \
                                                "edep_bending_screen": np.float64, "edep_bending_screen2": np.float64, \
                                                "edep_C1": np.float64, "edep_C2": np.float64, "edep_C3": np.float64, \
                                                "edep_C4": np.float64, "edep_C5": np.float64, "edep_C6": np.float64, \
                                                "edep_C7": np.float64, "edep_C8": np.float64, "edep_C9": np.float64, \
                                                "edep_C10": np.float64, "edep_C11": np.float64, "edep_C12": np.float64, \
                                                "edep_C13": np.float64, "edep_C14": np.float64, "edep_C15": np.float64, \
                                                "edep_C16": np.float64, "edep_C17": np.float64, "edep_C18": np.float64, \
                                                "edep_C19": np.float64, "edep_C20": np.float64, "edep_C21": np.float64, \
                                                "edep_C22": np.float64, "edep_C23": np.float64, "edep_C24": np.float64, \
                                                "edep_C25": np.float64
                                               })
        tree_out.extend({
                         "eventID": df_out_merged.eventID.values, \
                         "Tracker_NHit_X_1": df_out_merged.Tracker_NHit_X_1.values, "Tracker_NHit_Y_1": df_out_merged.Tracker_NHit_Y_1.values, \
                         "Tracker_NHit_X_2": df_out_merged.Tracker_NHit_X_2.values, "Tracker_NHit_Y_2": df_out_merged.Tracker_NHit_Y_2.values, \
                         "Tracker_X_1": df_out_merged.Tracker_X_1.values, "Tracker_Y_1": df_out_merged.Tracker_Y_1.values, \
                         "Tracker_X_2": df_out_merged.Tracker_X_2.values, "Tracker_Y_2": df_out_merged.Tracker_Y_2.values, \
                         "Ekin": df_out_merged.Ekin.values, \
                         "edep_APC1": df_out_merged.edep_APC1.values, "edep_APC2": df_out_merged.edep_APC2.values, \
                         "edep_calo": df_out_merged.edep_calo.values, "edep_screen": df_out_merged.edep_screen.values, \
                         "edep_calo2": df_out_merged.edep_calo2.values if 'edep_calo2' in df_out_merged.columns else np.nan , \
                         "edep_calo3": df_out_merged.edep_calo3.values if 'edep_calo3' in df_out_merged.columns else np.nan, \
                         "edep_bending_screen": df_out_merged.edep_bending_screen.values if 'edep_bending_screen' in df_out_merged.columns else np.nan, \
                         "edep_bending_screen2": df_out_merged.edep_bending_screen2.values if 'edep_bending_screen2' in df_out_merged.columns else np.nan, \
                         "edep_C1": df_out_merged.edep_C1.values if 'edep_C1' in df_out_merged.columns else np.nan, \
                         "edep_C2": df_out_merged.edep_C2.values if 'edep_C2' in df_out_merged.columns else np.nan, \
                         "edep_C3": df_out_merged.edep_C3.values if 'edep_C3' in df_out_merged.columns else np.nan, \
                         "edep_C4": df_out_merged.edep_C4.values if 'edep_C4' in df_out_merged.columns else np.nan, \
                         "edep_C5": df_out_merged.edep_C5.values if 'edep_C5' in df_out_merged.columns else np.nan, \
                         "edep_C6": df_out_merged.edep_C6.values if 'edep_C6' in df_out_merged.columns else np.nan, \
                         "edep_C7": df_out_merged.edep_C7.values if 'edep_C7' in df_out_merged.columns else np.nan, \
                         "edep_C8": df_out_merged.edep_C8.values if 'edep_C8' in df_out_merged.columns else np.nan, \
                         "edep_C9": df_out_merged.edep_C9.values if 'edep_C9' in df_out_merged.columns else np.nan, \
                         "edep_C10": df_out_merged.edep_C10.values if 'edep_C10' in df_out_merged.columns else np.nan, \
                         "edep_C11": df_out_merged.edep_C11.values if 'edep_C11' in df_out_merged.columns else np.nan, \
                         "edep_C12": df_out_merged.edep_C12.values if 'edep_C12' in df_out_merged.columns else np.nan, \
                         "edep_C13": df_out_merged.edep_C13.values if 'edep_C13' in df_out_merged.columns else np.nan, \
                         "edep_C14": df_out_merged.edep_C14.values if 'edep_C14' in df_out_merged.columns else np.nan, \
                         "edep_C15": df_out_merged.edep_C15.values if 'edep_C15' in df_out_merged.columns else np.nan, \
                         "edep_C16": df_out_merged.edep_C16.values if 'edep_C16' in df_out_merged.columns else np.nan, \
                         "edep_C17": df_out_merged.edep_C17.values if 'edep_C17' in df_out_merged.columns else np.nan, \
                         "edep_C18": df_out_merged.edep_C18.values if 'edep_C18' in df_out_merged.columns else np.nan, \
                         "edep_C19": df_out_merged.edep_C19.values if 'edep_C19' in df_out_merged.columns else np.nan, \
                         "edep_C20": df_out_merged.edep_C20.values if 'edep_C20' in df_out_merged.columns else np.nan, \
                         "edep_C21": df_out_merged.edep_C21.values if 'edep_C21' in df_out_merged.columns else np.nan, \
                         "edep_C22": df_out_merged.edep_C22.values if 'edep_C22' in df_out_merged.columns else np.nan, \
                         "edep_C23": df_out_merged.edep_C23.values if 'edep_C23' in df_out_merged.columns else np.nan, \
                         "edep_C24": df_out_merged.edep_C24.values if 'edep_C24' in df_out_merged.columns else np.nan, \
                         "edep_C25": df_out_merged.edep_C25.values if 'edep_C25' in df_out_merged.columns else np.nan, \
                        })
        tree_scr = rf_merged.mktree("scoringScreen", {
                                                "eventID": np.int32, "particle": np.int32, \
                                                "x": np.float64, "y": np.float64, "z": np.float64, \
                                                "px": np.float64, "py": np.float64, "pz": np.float64, \
                                                "t": np.float64, "E": np.float64, "parentID": np.int32, \
                                                "trackID": np.int32, "detID": np.int32
                                               })
        # conversion of string columns into integer columns (without this an error occurs)
        df_scr_merged["particle"] = df_scr_merged["particle"].astype(str) 
        part_dict = {item: i for i,item in enumerate(df_scr_merged.particle.unique())}
        df_scr_merged["particle"] = df_scr_merged["particle"].replace(part_dict)    
        tree_scr.extend({
                         "eventID": df_scr_merged.eventID.values, "particle": df_scr_merged.particle.values, \
                         "x": df_scr_merged.x.values, "y": df_scr_merged.y.values, "z": df_scr_merged.z.values, \
                         "px": df_scr_merged.px.values, "py": df_scr_merged.py.values, "pz": df_scr_merged.pz.values, \
                         "t": df_scr_merged.t.values, "E": df_scr_merged.E.values, "parentID": df_scr_merged.parentID.values, \
                         "trackID": df_scr_merged.trackID.values if 'trackID' in df_scr.columns else np.nan, \
                         "detID": df_scr_merged.detID.values if 'detID' in df_scr.columns else np.nan
                       })
    # print summary and return merged dataframes
    print("\n")
    if save_result:
        print("part_dict merged:", part_dict)
    else:
        part_dict = {}
    print("events_read:", events_read)
    print('%d files merged!\n' % (nfiles))
    return df_out_merged, df_scr_merged, part_dict


def merge_TestBeamOC_files(data_path, beVerbose=False, save_result=False):
    """
    Optimized (DeepSeek) function to merge root files of Geant4 simulations.
    """
    import uproot
    import numpy as np
    import pandas as pd
    import os
    
    # Filter for .root files only
    file_list = sorted([f for f in os.listdir(data_path) if f.endswith('.root')])
    if not file_list:
        print("No .root files found in directory")
        return pd.DataFrame(), pd.DataFrame()
    
    if beVerbose:
        print(f"Number of files to merge: {len(file_list)}")
    
    # Precompute file paths
    root_files = [os.path.join(data_path, f) for f in file_list]
    
    # Initialize lists for dataframes
    dataframes_out_list = []
    dataframes_scr_list = []
    
    events_read = 0
    nfiles = 0
    
    # Pre-defined column lists for faster checking
    expected_columns = [
        'edep_calo2', 'edep_calo3', 'edep_bending_screen', 'edep_bending_screen2',
        'edep_C1', 'edep_C2', 'edep_C3', 'edep_C4', 'edep_C5',
        'edep_C6', 'edep_C7', 'edep_C8', 'edep_C9', 'edep_C10',
        'edep_C11', 'edep_C12', 'edep_C13', 'edep_C14', 'edep_C15',
        'edep_C16', 'edep_C17', 'edep_C18', 'edep_C19', 'edep_C20',
        'edep_C21', 'edep_C22', 'edep_C23', 'edep_C24', 'edep_C25'
    ]
    
    for i, (root_file, filename) in enumerate(zip(root_files, file_list)):
        if beVerbose:
            print(f'Processing {filename} ({i+1}/{len(file_list)})...')
        
        try:
            # Open ROOT file
            rf = uproot.open(root_file)
            
            # Load both dataframes in parallel (faster for multiple trees)
            df_out = rf['outData'].arrays(library='pd')
            df_scr = rf['scoringScreen'].arrays(library='pd')
            
            # Get number of events (simpler method)
            n_events = df_scr['eventID'].nunique()
            
            # Vectorized eventID correction (much faster than per-event loop)
            event_offset = events_read
            df_out['eventID'] += event_offset
            df_scr['eventID'] += event_offset
            
            # Store dataframes
            dataframes_out_list.append(df_out)
            dataframes_scr_list.append(df_scr)
            
            # Update counters
            events_read += n_events
            nfiles += 1
            
            if beVerbose:
                print(f"  Events: {n_events}, Total: {events_read}")
                
        except Exception as e:
            print(f"Error processing {filename}: {e}")
            continue
    
    # Early return if no files processed
    if nfiles == 0:
        print("No files were successfully processed")
        return pd.DataFrame(), pd.DataFrame()
    
    # Merge dataframes (single concat operation)
    df_out_merged = pd.concat(dataframes_out_list, ignore_index=True, copy=False)
    df_scr_merged = pd.concat(dataframes_scr_list, ignore_index=True, copy=False)
    
    # Save results if requested
    if save_result and nfiles > 0:
        # Generate output filename efficiently
        last_file = file_list[-1].replace(".root", "")
        if '_' in last_file:
            # Find last underscore more efficiently
            parts = last_file.split('_')
            if len(parts) > 1:
                outputfile = '_'.join(parts[:-1]) + f'_merged{nfiles}files'
            else:
                outputfile = last_file + f'_merged{nfiles}files'
        else:
            outputfile = f'merged{nfiles}files'
        
        # Create output directory path
        output_dir = os.path.dirname(data_path)
        output_root = os.path.join(output_dir, '../', outputfile + '.root')
        
        # Create merged ROOT file
        rf_merged = uproot.recreate(output_root)
        
        # Prepare outData tree with optimized column handling
        # Get all possible columns from merged dataframe
        out_cols_present = [col for col in df_out_merged.columns if col in expected_columns]
        
        # Create dictionary for tree extension
        extend_dict = {
            "eventID": df_out_merged.eventID.values,
            "Tracker_NHit_X_1": df_out_merged.Tracker_NHit_X_1.values,
            "Tracker_NHit_Y_1": df_out_merged.Tracker_NHit_Y_1.values,
            "Tracker_NHit_X_2": df_out_merged.Tracker_NHit_X_2.values,
            "Tracker_NHit_Y_2": df_out_merged.Tracker_NHit_Y_2.values,
            "Tracker_X_1": df_out_merged.Tracker_X_1.values,
            "Tracker_Y_1": df_out_merged.Tracker_Y_1.values,
            "Tracker_X_2": df_out_merged.Tracker_X_2.values,
            "Tracker_Y_2": df_out_merged.Tracker_Y_2.values,
            "Ekin": df_out_merged.Ekin.values,
            "edep_APC1": df_out_merged.edep_APC1.values,
            "edep_APC2": df_out_merged.edep_APC2.values,
            "edep_calo": df_out_merged.edep_calo.values,
            "edep_screen": df_out_merged.edep_screen.values,
        }
        
        # Add optional columns efficiently
        for col in out_cols_present:
            extend_dict[col] = df_out_merged[col].values
        
        # Create tree definition dynamically
        tree_def = {
            "eventID": np.int32,
            "Tracker_NHit_X_1": np.int32, "Tracker_NHit_Y_1": np.int32,
            "Tracker_NHit_X_2": np.int32, "Tracker_NHit_Y_2": np.int32,
            "Tracker_X_1": np.float64, "Tracker_Y_1": np.float64,
            "Tracker_X_2": np.float64, "Tracker_Y_2": np.float64,
            "Ekin": np.float64,
            "edep_APC1": np.float64, "edep_APC2": np.float64,
            "edep_calo": np.float64, "edep_screen": np.float64,
        }
        
        # Add optional columns to tree definition
        for col in out_cols_present:
            tree_def[col] = np.float64
        
        tree_out = rf_merged.mktree("outData", tree_def)
        tree_out.extend(extend_dict)
        
        # Prepare scoringScreen tree
        # Convert particle strings to categorical codes efficiently
        part_dict = {item: i for i,item in enumerate(df_scr_merged.particle.unique())}
        df_scr_merged["particle"] = pd.Categorical(df_scr_merged["particle"]).codes
        
        # Create extend dictionary for scoringScreen
        scr_extend_dict = {
            "eventID": df_scr_merged.eventID.values,
            "particle": df_scr_merged.particle.values,
            "x": df_scr_merged.x.values,
            "y": df_scr_merged.y.values,
            "z": df_scr_merged.z.values,
            "px": df_scr_merged.px.values,
            "py": df_scr_merged.py.values,
            "pz": df_scr_merged.pz.values,
            "t": df_scr_merged.t.values,
            "E": df_scr_merged.E.values,
            "parentID": df_scr_merged.parentID.values,
        }
        
        # Add optional columns if they exist
        if 'trackID' in df_scr_merged.columns:
            scr_extend_dict["trackID"] = df_scr_merged.trackID.values
        else:
            scr_extend_dict["trackID"] = np.full(len(df_scr_merged), -1, dtype=np.int32)
            
        if 'detID' in df_scr_merged.columns:
            scr_extend_dict["detID"] = df_scr_merged.detID.values
        else:
            scr_extend_dict["detID"] = np.full(len(df_scr_merged), -1, dtype=np.int32)
        
        # Define tree structure
        tree_scr_def = {
            "eventID": np.int32, "particle": np.int32,
            "x": np.float64, "y": np.float64, "z": np.float64,
            "px": np.float64, "py": np.float64, "pz": np.float64,
            "t": np.float64, "E": np.float64, "parentID": np.int32,
            "trackID": np.int32, "detID": np.int32
        }
        
        tree_scr = rf_merged.mktree("scoringScreen", tree_scr_def)
        tree_scr.extend(scr_extend_dict)
            
    # Print summary and return merged dataframes
    print("\n")
    if save_result:
        print("part_dict merged:", part_dict)
    else:
        part_dict = {}
    print(f"Total events read: {events_read}")
    print(f"Files successfully merged: {nfiles}\n")
    
    return df_out_merged, df_scr_merged, part_dict
#######################################################################################################
