# Set of functions to read (and plot) Edep in various volumes, scored in different ways.
# These functions were conceived for (oriented) calorimeters, thus Edep is scored in GeV.  

# Import the required libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#%matplotlib ipympl
import os
import uproot
import math

from G4_utils import * #my custom functions of general utility



def get_scored_Edep_from_txt(path, set_list, add_to_all_pre, add_to_all_post, normEvents):
    """
    Function to read the total Edep (saved in txt files) for many runs (given by set_list).
    Pass always a list (set_list), even with if you have only one element.
    For each run, it foresees a file with 2 header rows (the 1st with Nevents and the 2nd
    with the column labels) and and many columns with volumeID, elementID, Edep, stdEdep.   
    It returns [volumeID, elementID, Edep, stdEdep], namely a list with each of element 
    being a dictionary with items of set_list as keys and values suggested by the name.
    """
    
    volumeID = {}
    elementID = {}
    Edep = {}
    stdEdep = {}
    
    for item in set_list:
        filename = path + add_to_all_pre + item + add_to_all_post
        
        with open(filename) as f:
            all_data = [line.strip() for line in f.readlines()]
            Nevents = int(all_data[0].split(':')[1])
            columns = all_data[1].split(' ')
            data = all_data[2:]
            N = len(data)
            
            if normEvents:
                cf = 1. / Nevents
            else:
                cf = 1.

            volumeID_list = []
            elementID_list = []
            Edep_list = []
            stdEdep_list = []
            for k in range(N):
                temp = data[k].split(' ')
                volumeID_list.append(int(temp[0]) )
                elementID_list.append(int(temp[1]))
                Edep_list.append(float(temp[2]) * cf)
                stdEdep_list.append(float(temp[3]) * np.sqrt(cf))
            
            volumeID[item] = volumeID_list
            elementID[item] = elementID_list
            Edep[item] = Edep_list
            stdEdep[item] = stdEdep_list
            
    return volumeID, elementID, Edep, stdEdep



def read_and_plot_Edep_thread(filename, beVerbose=True, doPlots=True, saveFigs=False):
    """
    Function to read (and plot) Edep accumulated (through built-in scorer) 
    per thread in two volumes (PS and CL).
    It returns [df_PS_thread_edep, df_CL_thread_edep, edep_sum_PS, edep_sum_CL].
    It does plots [by default].
    """
    
    path_split = filename.split('/')[1:-1]
    path = ""
    for item in path_split:
        path = path + '/' + item
    path = path + '/' 
    
    rf = uproot.open(filename)
    branches = ["crystalID", "edep"]
    df_PS_thread = rf['PreShowerCrystals_Edep_thread'].arrays(branches, library='pd')
    df_CL_thread = rf['CalorimeterCrystals_Edep_thread'].arrays(branches, library='pd') 
    
    crystalID_PS = df_PS_thread.crystalID.unique()
    crystalID_CL = df_CL_thread.crystalID.unique()
    
    Ncrystals_PS = len(crystalID_PS)
    Ncrystals_CL = len(crystalID_CL)
    Nthreads = int(len(df_CL_thread)/len(crystalID_CL))
    
    edep_PS = {}
    for ID in crystalID_PS:
        edep = df_PS_thread.loc[df_PS_thread['crystalID'] == ID].edep.values
        edep_PS[ID] = edep
    edep_CL = {}
    for ID in crystalID_CL:
        edep = df_CL_thread.loc[df_CL_thread['crystalID'] == ID].edep.values
        edep_CL[ID] = edep
    
    df_PS_thread_edep = pd.DataFrame(edep_PS)
    df_CL_thread_edep = pd.DataFrame(edep_CL)
    
    df_PS_edep_stat = df_PS_thread_edep.describe()
    df_CL_edep_stat = df_CL_thread_edep.describe()
    
    edep_mean_PS = []
    edep_std_PS = []
    edep_sum_PS = []
    if beVerbose:
        print("Energy deposited per thread in PreShower crystals:")
    for ID in crystalID_PS:
        edep_mean_PS.append(df_PS_edep_stat[ID].values[1])
        edep_std_PS.append(df_PS_edep_stat[ID].values[2])
        edep_sum_PS.append(sum(edep_PS[ID]))
        if beVerbose:
            print("Edep[%d] = %.4f +/- %.4f GeV" % (ID, edep_mean_PS[ID], edep_std_PS[ID]))
    print("\n")
    edep_mean_CL = []
    edep_std_CL = []
    edep_sum_CL = []
    if beVerbose:
        print("Energy deposited per thread in Calorimeter crystals:")
    for ID in crystalID_CL:
        edep_mean_CL.append(df_CL_edep_stat[ID].values[1])
        edep_std_CL.append(df_CL_edep_stat[ID].values[2])
        edep_sum_CL.append(sum(edep_CL[ID]))
        if beVerbose:
            print("Edep[%d] = %.4f +/- %.4f GeV" % (ID, edep_mean_CL[ID], edep_std_CL[ID]))
    print("\n")
    
    total_Edep_PS = sum(edep_sum_PS)
    print('total_Edep_PS: %.0f GeV' % total_Edep_PS)
    total_Edep_CL = sum(edep_sum_CL)
    print('total_Edep_CL: %.0f GeV\n\n' % total_Edep_CL)
    
    if doPlots:
        fig = plt.figure(figsize=(6, 6))
        fs = 14
        ms = 8
        plt.subplot(2,1,1)
        edep_image_PS = np.expand_dims(edep_sum_PS/total_Edep_PS, axis=0)
        plt.imshow(10*np.log10(edep_image_PS), cmap='plasma')
        cbar = plt.colorbar(location="bottom")
        cbar.set_label('Edep (dB)', fontsize=fs, rotation=0)
        cbar.ax.tick_params(labelsize=fs)
        plt.title('PreShower crystals', fontsize=fs)
        plt.xlabel("", fontsize=fs)
        plt.ylabel("", fontsize=fs)
        plt.xticks(fontsize=fs, rotation=0)
        plt.yticks([],fontsize=fs, rotation=0)
        plt.subplot(2,1,2)
        edep_image_CL = np.expand_dims(edep_sum_CL/total_Edep_CL, axis=0)
        plt.imshow(10*np.log10(edep_image_CL), cmap='plasma')
        cbar = plt.colorbar(location="bottom")
        cbar.set_label('Edep (dB)', fontsize=fs, rotation=0)
        cbar.ax.tick_params(labelsize=fs)
        plt.title('Calorimeter crystals', fontsize=fs)
        plt.xlabel("", fontsize=fs)
        plt.ylabel("", fontsize=fs)
        plt.xticks(fontsize=fs, rotation=0)
        plt.yticks([],fontsize=fs, rotation=0)   
        if saveFigs:
            plt.savefig(path + 'Edep_crystals.jpg')
        plt.tight_layout(pad=1)
        plt.show()  
        #plt.close()
    
    return df_PS_thread_edep, df_CL_thread_edep, edep_sum_PS, edep_sum_CL



def read_and_plot_Edep_event(filename, normEvents, beVerbose=True, doPlots=True, saveFigs=False):
    """
    Function to read (and plot) Edep accumulated (through built-in scorer) 
    per event in two volumes (PS and CL).
    It returns [edep_PS_event, edep_CL_event, edep_PS, edep_CL, df_PS_event, df_CL_event].
    It does plots [by default].
    """
    
    path_split = filename.split('/')[1:-1]
    path = ""
    for item in path_split:
        path = path + '/' + item
    path = path + '/'  
    
    rf = uproot.open(filename)
    branches = ["eventID", "crystalID", "edep"]
    df_PS_event = rf['EdepPreShowerCrystals'].arrays(branches, library='pd')
    df_CL_event = rf['EdepCalorimeterCrystals'].arrays(branches, library='pd') 
    
    crystalID_PS = df_PS_event.crystalID.unique()
    crystalID_CL = df_CL_event.crystalID.unique()

    Ncrystals_PS = len(crystalID_PS)
    Ncrystals_CL = len(crystalID_CL)
    Nevents = int(len(df_CL_event)/len(crystalID_CL))
    if normEvents:
        nf = Nevents
        if beVerbose:
            print('ATTENTION: Edep is intended per event')
    else:
        nf = 1

    edep_PS_event = {}
    for ID in crystalID_PS:
        edep = df_PS_event.loc[df_PS_event['crystalID'] == ID].edep.values
        edep_PS_event[ID] = edep
    edep_CL_event = {}
    for ID in crystalID_CL:
        edep = df_CL_event.loc[df_CL_event['crystalID'] == ID].edep.values
        edep_CL_event[ID] = edep

    edep_PS = []
    total_edep_PS = 0
    for ID in edep_PS_event.keys():
        edep_PS.append(np.sum(edep_PS_event[ID])/nf)
        if beVerbose:
            print('PS_Edep[%d]: %.6f GeV' % (ID, edep_PS[ID]))
        total_edep_PS = total_edep_PS + edep_PS[ID]
    edep_CL = []
    total_edep_CL = 0
    for ID in edep_CL_event.keys():
        edep_CL.append(np.sum(edep_CL_event[ID])/nf)
        if beVerbose:
            print('CL_Edep[%d]: %.6f GeV' % (ID, edep_CL[ID]))
        total_edep_CL = total_edep_CL + edep_CL[ID]
        
    if beVerbose:
        print('total_Edep in PreSHower: %.2f GeV' % (total_edep_PS))
        print('total_Edep in Calorimeter: %.2f GeV' % (total_edep_CL))
        print('\n')

    if doPlots:
        edep_image_PS = np.expand_dims(edep_PS/total_edep_PS, axis=0)
        if "layers" in filename:
            CL_layers = int(filename.split("layers")[0][-1])
            print("CL_layers:", CL_layers)
            edep_image_CL = np.reshape(edep_CL/total_edep_CL, 
                                       (CL_layers, int(Ncrystals_CL/CL_layers)))
            ylbl = "layers"
            yticks = list(np.linspace(0, CL_layers-1, CL_layers, dtype=int))
            yticks_lbl = [str(item) for item in list(np.linspace(1, CL_layers, CL_layers, dtype=int))]
        else:
            edep_image_CL = np.expand_dims(edep_CL/total_edep_CL, axis=0)
            ylbl = ""
            yticks = []
            yticks_lbl = []
        
        fig = plt.figure(figsize=(6, 6))
        fs = 14
        ms = 8
        plt.subplot(2,1,1)
        plt.imshow(10*np.log10(edep_image_PS), cmap='plasma')
        cbar = plt.colorbar(location="bottom")
        cbar.set_label('Edep (dB)', fontsize=fs, rotation=0)
        cbar.ax.tick_params(labelsize=fs)
        plt.title('PreShower crystals', fontsize=fs)
        plt.xlabel("", fontsize=fs)
        plt.ylabel("", fontsize=fs)
        plt.xticks(fontsize=fs, rotation=0)
        plt.yticks([],fontsize=fs, rotation=0)
        plt.subplot(2,1,2)
        plt.imshow(10*np.log10(edep_image_CL), cmap='plasma')
        cbar = plt.colorbar(location="bottom")
        cbar.set_label('Edep (dB)', fontsize=fs, rotation=0)
        cbar.ax.tick_params(labelsize=fs)
        plt.title('Calorimeter crystals', fontsize=fs)
        plt.xlabel("", fontsize=fs)
        plt.ylabel(ylbl, fontsize=fs)
        plt.xticks(fontsize=fs, rotation=0)
        plt.yticks(ticks=yticks, labels=yticks_lbl, fontsize=fs, rotation=0)         
        if saveFigs:
            plt.savefig(path + 'Edep_crystals.jpg')
        plt.tight_layout(pad=1)
        plt.show()  
        #plt.close()
    
    return edep_PS_event, edep_CL_event, edep_PS, edep_CL, df_PS_event, df_CL_event



def read_Edep_event(filename, normEvents, beVerbose=True):
    """
    Function to read Edep accumulated (through SteppingAction) 
    per event in various Volumes. It is general.
    It returns [edep_event, edep, df_event, Nvolumes, Nevents].
    """
    
    rf = uproot.open(filename)
    df_event = rf['Edep'].arrays(library='pd')

    volume_IDs = df_event.volumeID.unique()
    
    Nvolumes = len(volume_IDs)
    Nevents = int(len(df_event)/Nvolumes)
    #print("Nevents:", Nevents)
    if normEvents:
        nf = Nevents
        if beVerbose:
            print('ATTENTION: Edep is intended per event')
    else:
        nf = 1
        
    edep_event = {}
    for ID in volume_IDs:
        edep = df_event.loc[df_event['volumeID'] == ID].edep.values
        edep_event[ID] = edep
        
    edep = []
    total_edep = 0
    for ID in edep_event.keys():
        edep.append(np.sum(edep_event[ID])/nf)
        if beVerbose:
            print('Edep[%d]: %.6f GeV' % (ID, edep[ID]))
        total_edep = total_edep + edep[ID]
    if beVerbose:
        print('total_Edep in Volumes: %.2f GeV' % (total_edep))
        print('\n')
        
    return edep_event, edep, df_event, Nvolumes, Nevents



def read_and_plot_Edep_VoxelScorer(filename, normEvents, Nevents, 
                                   tranvsizeX, tranvsizeY, tranvsizeZ, \
                                   axis_proj, proj_type, IWantPlotVoxel=True, \
                                   p2m=3, beVerbose=True, doPlots=True, \
                                   fs=26, cm='viridis', saveFigs=False):
    """
    Function to read (and plot) the total Edep accumulated through 
    a custom VoxelScorer associated to a volume. It is general.
    It returns [SLICES, proj, profileZ, (x,y,z), profileR, cum_profileR, r_list].
    """
    
    # infere the path of the file (for figure saving)
    path_split = filename.split('/')[1:-1]
    path = ""
    for item in path_split:
        path = path + '/' + item
    path = path + '/'

    # read the file
    temp_list = []
    full_list = []
    with open(filename, 'r') as f:
        header = f.readline()
        nx, ny, nz = [int(item) for item in header.split('\n')[0].split(' ')]
        for line in f:
            temp_list = []       
            temp = line.split(' \n')[0].split(' ')
            for item in temp:
                temp_list.append(float(item))
            full_list.append(temp_list)
    full_array = np.array(full_list)

    # print the total value
    total = np.sum(full_array)
    if beVerbose:
        print('total Edep (before normalization): %.2f GeV' % total)

    # set title string
    if "PreShower" in filename:
        ttlstr = 'PreShower'
    elif "Calorimeter" in filename:
        ttlstr = 'Calorimeter'
    else:
        ttlstr = ''

    # normalize edep per event
    if normEvents:
        full_array = full_array / Nevents
        clblabel='Edep per event (GeV)'
    else:
        clblabel='Edep (GeV)'

    # generate 3D data matrix
    SLICES = full_array.reshape(nz, nx, ny)
    
    # set geometrical size
    Xmin = -tranvsizeX*0.5        
    Xmax = tranvsizeX*0.5         
    Ymin = -tranvsizeY*0.5        
    Ymax = tranvsizeY*0.5
    Zmin = 0
    Zmax = tranvsizeZ

    # generate vectors of points linearly spaced (corrected in 24/12/2023)
    Xedges = np.linspace(Xmin, Xmax, nx+1) #[mm]
    Yedges = np.linspace(Ymin, Ymax, ny+1) #[mm]
    Zedges = np.linspace(Zmin, Zmax, nz+1) #[mm]
    dx = Xedges[1] - Xedges[0]
    dy = Yedges[1] - Yedges[0]
    dz = Zedges[1] - Zedges[0]
    #print("dx:", dx, "mm, dy:", dy, "mm, dz:", dz, "mm")
    x = Xedges[:-1] + dx*0.5
    y = Yedges[:-1] + dy*0.5
    z = Zedges[:-1] + dz*0.5  

    # mean/sum distribution in the desired projection
    if proj_type == 1:
        proj = np.mean(SLICES, axis=axis_proj)
        projstr = 'projection type: mean'
    else:
        proj = np.sum(SLICES, axis=axis_proj)
        projstr = 'projection type: sum'
    
    # set x,y labels
    if axis_proj == 0:
        if beVerbose:
            print('projection: axial')
        if IWantPlotVoxel:
            xlbl='X (voxel)'
            ylbl='Y (voxel)'
        else:
            xlbl='X (mm)'
            ylbl='Y (mm)'
        X, Y = np.meshgrid(x, y)
        proj = np.transpose(proj)
    elif axis_proj == 1:
        if beVerbose:
            print('projection: sagittal')  
        if IWantPlotVoxel:
            xlbl='Y (voxel)'
            ylbl='Z (voxel)'
        else:            
            xlbl='Y (mm)'
            ylbl='Z (mm)'
        X, Y = np.meshgrid(y, z)
    else:
        if beVerbose:
            print('projection: coronal') 
        if IWantPlotVoxel:
            xlbl='X (voxel)'
            ylbl='Z (voxel)'
        else:
            xlbl='X (mm)'
            ylbl='Z (mm)'
        X, Y = np.meshgrid(x, z)

    # plot the desired projection
    if (doPlots and X.shape[0]>1 and X.shape[1]>1):
        print(projstr);
        print('projection shape:', proj.shape, 'voxel')      
        
        fig = plt.figure(figsize=(10, 8))
        Nlev = 50
        if IWantPlotVoxel:
            #plt.imshow(proj, origin='upper', cmap=cm)
            plt.contourf(proj, Nlev, cmap=cm)
            #plt.gca().set_aspect('equal')
        else:
            plt.contourf(X, Y, proj, Nlev, cmap=cm)
        plt.gca().invert_yaxis()
        cbar = plt.colorbar()
        cbar.ax.tick_params(labelsize=fs)
        cbar.set_label(clblabel, fontsize=fs, rotation=90)
        plt.title(ttlstr, fontsize=fs)
        plt.xlabel(xlbl, fontsize=fs)
        plt.ylabel(ylbl, fontsize=fs)
        plt.xticks(fontsize=fs, rotation=0)
        plt.yticks(fontsize=fs, rotation=0)
        if saveFigs:
            plt.savefig(path + 'Edep_proj_voxelscorer.jpg')
        plt.show()  
        #plt.close()

    # calculate the profile along Z direction
    profileZ = np.zeros(proj.shape[0])
    if axis_proj > 0:
        proj_4profZ = proj
    else:
        print('Since axial projection was selected, profileZ is arbitrarly calculted w.r.t coronal projection.')
        if proj_type == 1:
            proj_4profZ = np.mean(SLICES, axis=2)
        else:
            proj_4profZ = np.sum(SLICES, axis=2)
    ic = int(proj.shape[1]/2+1)
    if proj_type == 1:
        for j in range(nz):
            profileZ[j] = np.mean(proj_4profZ[j][ic-p2m:ic+p2m])
            #profileZ[j] = np.mean(proj_4profZ[j])
    else:
        for j in range(nz):
            profileZ[j] = np.sum(proj_4profZ[j][ic-p2m:ic+p2m])
            #profileZ[j] = np.sum(proj_4profZ[j])
                
    # calculate the radial profile
    if proj_type == 1:
        axial_proj = np.mean(SLICES, axis=0)
    else:
        axial_proj = np.sum(SLICES, axis=0)
    axial_proj = np.transpose(axial_proj)
    XY = np.meshgrid(x, y)
    R = np.sqrt(XY[0]*XY[0] + XY[1]*XY[1])
    r_list = np.linspace(0.001, np.max(R), 200)
    profileR = np.zeros(r_list.shape)
    cum_profileR = np.zeros(r_list.shape)
    for i,r in enumerate(r_list):
        sel = list(axial_proj[R < r])
        profileR[i] = np.nanmean(sel)
        cum_profileR[i] = np.sum(sel)
        
    # return a set of variables
    return SLICES, proj, profileZ, (x,y,z), profileR, cum_profileR, r_list



def read_VoxelScorer(filename, normEvents, Nevents, \
                     tranvsizeX, tranvsizeY, tranvsizeZ, \
                     axis_proj, proj_type, \
                     p2m=3, beVerbose=True):
    """
    Function to read the quantity (Edep, Dose, Kerma) accumulated through 
    a custom VoxelScorer associated to a volume. It is general.
    It returns [SLICES, proj, profileZ, (x,y,z), profileR, cum_profileR, r_list].
    """
    
    # infere the path of the file (for figure saving)
    path_split = filename.split('/')[1:-1]
    path = ""
    for item in path_split:
        path = path + '/' + item
    path = path + '/'

    # read the file
    temp_list = []
    full_list = []
    with open(filename, 'r') as f:
        header = f.readline()
        nx, ny, nz = [int(item) for item in header.split('\n')[0].split(' ')]
        for line in f:
            temp_list = []       
            temp = line.split(' \n')[0].split(' ')
            for item in temp:
                temp_list.append(float(item))
            full_list.append(temp_list)
    full_array = np.array(full_list)

    # print the total value
    total = np.sum(full_array)
    if beVerbose:
        print('total value (before normalization): %.4e' % total)

    # normalize quantity per event
    if normEvents:
        full_array = full_array / Nevents

    # generate 3D data matrix
    SLICES = full_array.reshape(nz, nx, ny)
    
    # set geometrical size
    Xmin = -tranvsizeX*0.5        
    Xmax = tranvsizeX*0.5         
    Ymin = -tranvsizeY*0.5        
    Ymax = tranvsizeY*0.5
    Zmin = 0
    Zmax = tranvsizeZ

    # generate vectors of points linearly spaced (corrected in 24/12/2023)
    Xedges = np.linspace(Xmin, Xmax, nx+1) #[mm]
    Yedges = np.linspace(Ymin, Ymax, ny+1) #[mm]
    Zedges = np.linspace(Zmin, Zmax, nz+1) #[mm]
    dx = Xedges[1] - Xedges[0]
    dy = Yedges[1] - Yedges[0]
    dz = Zedges[1] - Zedges[0]
    #print("dx:", dx, "mm, dy:", dy, "mm, dz:", dz, "mm")
    x = Xedges[:-1] + dx*0.5
    y = Yedges[:-1] + dy*0.5
    z = Zedges[:-1] + dz*0.5  

    # mean/sum distribution in the desired projection
    if proj_type == 1:
        proj = np.mean(SLICES, axis=axis_proj)
        projstr = 'projection type: mean'
    else:
        proj = np.sum(SLICES, axis=axis_proj)
        projstr = 'projection type: sum'
    
    # set x,y labels
    if axis_proj == 0:
        if beVerbose:
            print('projection: axial')
        X, Y = np.meshgrid(x, y)
        proj = np.transpose(proj)
    elif axis_proj == 1:
        if beVerbose:
            print('projection: sagittal') 
        X, Y = np.meshgrid(y, z)
    else:
        if beVerbose:
            print('projection: coronal')
        X, Y = np.meshgrid(x, z)

    # print the projection shape
    print(projstr);
    print('projection shape:', proj.shape, 'voxel')      

    # calculate the profile along Z direction
    profileZ = np.zeros(proj.shape[0])
    if axis_proj > 0:
        proj_4profZ = proj
    else:
        print('Since axial projection was selected, profileZ is arbitrarly calculted w.r.t coronal projection.')
        if proj_type == 1:
            proj_4profZ = np.mean(SLICES, axis=2)
        else:
            proj_4profZ = np.sum(SLICES, axis=2)
    ic = int(proj.shape[1]/2+1)
    if proj_type == 1:
        for j in range(nz):
            profileZ[j] = np.mean(proj_4profZ[j][ic-p2m:ic+p2m])
            #profileZ[j] = np.mean(proj_4profZ[j])
    else:
        for j in range(nz):
            profileZ[j] = np.sum(proj_4profZ[j][ic-p2m:ic+p2m])
            #profileZ[j] = np.sum(proj_4profZ[j])
                
    # calculate the radial profile
    if proj_type == 1:
        axial_proj = np.mean(SLICES, axis=0)
    else:
        axial_proj = np.sum(SLICES, axis=0)
    axial_proj = np.transpose(axial_proj)
    XY = np.meshgrid(x, y)
    R = np.sqrt(XY[0]*XY[0] + XY[1]*XY[1])
    r_list = np.linspace(0.001, np.max(R), 200)
    profileR = np.zeros(r_list.shape)
    cum_profileR = np.zeros(r_list.shape)
    for i,r in enumerate(r_list):
        sel = list(axial_proj[R < r])
        profileR[i] = np.nanmean(sel)
        cum_profileR[i] = np.sum(sel)
        
    # return a set of variables
    return SLICES, proj, profileZ, (x,y,z), profileR, cum_profileR, r_list



def read_Edep_BoxMesh(filename, normEvents, Nevents, 
                      tranvsizeX, tranvsizeY, tranvsizeZ):
    """
    Function to read the total Edep accumulated through a built-in BoxMesh scorer.
    It is general. It returns [data, (x,y,z)], where:
    data is a dictionary with 10 columns defined as follows
    {"ind_x", "ind_y", "ind_z", "eDep", "x", "y", "z", "r" ,"eDep_err", "eDepDensity"}
    and (x,y,z) are the coordinates [mm] of the voxel centers.
    """
    
    def find(cond, N=1e7):    
        """
        function, equivalent to Matlab's one, that returns
        a list with the indices of an array/matrix
        that satisfy a given condition.
        """
        cond = cond.reshape(np.size(cond), 1).ravel()
        index = list(np.argwhere(cond).ravel())
        return index[:N]
    
    # retrieve data and put them into a dataframe
    data = pd.read_csv(
        filename, skiprows=3, names=["ind_x", "ind_y", "ind_z", "eDep", "eDep2", "Nentry"]
    )
    
    # calculate eDep uncertainty 
    data["eDep_err"] = (data["eDep2"]/data["Nentry"] - 
                        (data["eDep"]/data["Nentry"])**2)**0.5
    data.fillna(0, inplace=True)
    data = data.drop(['eDep2', 'Nentry'])

    # retrieve the number of voxels in each direction
    Nvoxel = len(data)
    zzzz = find(np.array(data["ind_z"]) == 0, 2)
    if len(zzzz) > 1:
        nz = zzzz[1]
    else:
        nz = len(data["ind_y"])    
    ny = np.max([round((np.max(find(np.array(data["ind_y"])==0, nz+1))-1)/nz), 1])
    nx = max(np.array(data["ind_x"]))+1
    
    # set geometrical size
    Xmin = -tranvsizeX*0.5        
    Xmax = tranvsizeX*0.5         
    Ymin = -tranvsizeY*0.5        
    Ymax = tranvsizeY*0.5
    Zmin = 0
    Zmax = tranvsizeZ

    # generate vectors of points linearly spaced (corrected in 24/12/2023)
    Xedges = np.linspace(Xmin, Xmax, nx+1) #[mm]
    Yedges = np.linspace(Ymin, Ymax, ny+1) #[mm]
    Zedges = np.linspace(Zmin, Zmax, nz+1) #[mm]
    dx = Xedges[1] - Xedges[0]
    dy = Yedges[1] - Yedges[0]
    dz = Zedges[1] - Zedges[0]
    #print("dx:", dx, "mm, dy:", dy, "mm, dz:", dz, "mm")
    x = Xedges[:-1] + dx*0.5
    y = Yedges[:-1] + dy*0.5
    z = Zedges[:-1] + dz*0.5

    # physical position inside the mesh data
    data["x"] = dx * (data["ind_x"] + 0.5) - tranvsizeX*0.5 #central x [mm]
    data["y"] = dy * (data["ind_y"] + 0.5) - tranvsizeY*0.5 #central y [mm]
    data["z"] = dz * (data["ind_z"] + 0.5)                  #central z [mm]
    data["r"] = np.sqrt(data["x"]**2 + data["y"]**2)        #central r [mm]
    
    # eDep Map
    if normEvents:
        data["eDep"] = data["eDep"] / Nevents
        data["eDep_err"] = data["eDep_err"] / Nevents
    data["eDepDensity"] = data["eDep"] / (dz * dy * dx) #[MeV/mm**3] or [MeV/(mm**3 event)]

    # return
    return data, (x,y,z)
