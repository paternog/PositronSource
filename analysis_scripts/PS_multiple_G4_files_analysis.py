#!/usr/bin/env python
# coding: utf-8


# Read the output of a set of Geant4 PositronSource simulations


## Import the required libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#%matplotlib ipympl
from matplotlib.ticker import AutoMinorLocator
import os
import uproot
import math
import re
from scipy import interpolate

from G4_utils import *
from G4_read_output_files import *

### required libraries to make plots according Mattia's style
from scipy.optimize import curve_fit
from matplotlib.colors import LogNorm
import succolib as sl # Mattia's package (install it with: pip install succolib)
### Set plot style ###
import mplhep as hep
import plot_settings as settings
style = settings.style #use this with "import plot_settings as settings"
style['figure.figsize'] = (style['figure.figsize'][0], 8.50) #manual tweak to figure height
plt.style.use(style)

# define MPLBACKEND
os.environ['MPLBACKEND'] = 'Agg' 
MPLBACKEND = os.environ['MPLBACKEND'] 
print('\nMPLBACKEND:', MPLBACKEND, '\n')

# tic
import time
start = time.time()


## Input and Settings
###############################################################################
## Input path and base filename
Ein = 2.86 #GeV
sigma_in = 1.0 #mm


#path = "/home/paterno/geant4-apps/PositronSource-build/output/conventional_statistical_analysis/"
#rad_th = list(np.linspace(0, 9, 10, dtype=int)) #used simply as iterator
#Dcm = [''] #cm
#conv_th = [17.6] #mm 

#path = "/home/paterno/geant4-apps/PositronSource-build/output/results_2.86GeV_sigma1.0mm_W5-20mm_conventional_gp/"
#rad_th = list(np.linspace(5.0, 20.0, 16)) #mm
#Dcm = [''] #cm
#conv_th = [''] #mm

#path = "/home/paterno/geant4-apps/PositronSource-build/output/results_2.86GeV_sigma1.0mm_W8-16mm_crystalline_gp/"
#rad_th = list(np.linspace(8.0, 16.0, 9)) #mm
#Dcm = [''] #cm
#conv_th = [''] #mm

#path = "C:/DATI/Geant4/PositronSource-build/output/results_2.86GeV_sigma1.0mm_W12mm_crystalline_mis_HT_gp/"
path = "/home/paterno/geant4-apps/PositronSource-build/output/results_2.86GeV_sigma1.0mm_W12mm_crystalline_mis_HT_gp/"
rad_th = ['12.0_0.000_0.065', '12.0_0.001_0.065', '12.0_0.002_0.065', '12.0_0.003_0.065', \
          '12.0_0.004_0.065', '12.0_0.005_0.065', '12.0_0.006_0.065', '12.0_0.007_0.065', \
          '12.0_0.008_0.065', '12.0_0.009_0.065', '12.0_0.010_0.065'
         ] #mm_rad_A
Dcm = [''] #cm
conv_th = [''] #mm

#path = "/home/paterno/geant4-apps/PositronSource-build/output/results_6GeV_sigma0.5mm_W4-24mm_conventional_gp/"
#rad_th = list(np.linspace(4.0, 24.0, 21)) #mm
#Dcm = [''] #cm
#conv_th = [''] #mm

#path = "/home/paterno/geant4-apps/PositronSource-build/output/results_6GeV_sigma0.5mm_W9-15mm_crystalline_gp/"
#rad_th = list(np.linspace(9.0, 15.0, 7)) #mm
#Dcm = [''] #cm
#conv_th = [''] #mm

#path = "/home/paterno/geant4-apps/PositronSource-build/output/results_6GeV_sigma0.5mm_W12mm_crystalline_mis_HT_gp/"
#rad_th = ['12.0_0.000_0.065', '12.0_0.001_0.065', '12.0_0.002_0.065', '12.0_0.003_0.065', \
#          '12.0_0.004_0.065', '12.0_0.005_0.065', '12.0_0.006_0.065', '12.0_0.007_0.065', \
#          '12.0_0.008_0.065', '12.0_0.009_0.065', '12.0_0.010_0.065'
#         ] #mm_rad_A
#Dcm = [''] #cm
#conv_th = [''] #mm

#path = "/home/paterno/geant4-apps/PositronSource-build/output/results_6GeV_sigma0.5mm_W_radiator1-2mm_D0cm_target11.6mm_gp/"
#rad_th = list(np.linspace(1.0, 2.0, 11)) #mm
#Dcm = [0] #cm
#conv_th = [11.6] #mm

#path = "/home/paterno/geant4-apps/PositronSource-build/output/results_6GeV_sigma0.5mm_W_radiator2mm_D0-60cm_target9mm_gp/"
#rad_th = [2.0] #mm
#Dcm = list(np.linspace(0, 60, 13, dtype=int)) #cm
#conv_th = [9.0] #mm

#path = "/home/paterno/geant4-apps/PositronSource-build/output/results_6GeV_sigma0.5mm_W_radiator2mm_D0cm_target6-12mm_gp/"
#rad_th = [2.0] #mm
#Dcm = [0] #cm
#conv_th = list(np.linspace(6.0, 12.0, 7)) #mm

#path = "/home/paterno/geant4-apps/PositronSource-build/output/results_6GeV_sigma0.5mm_W_radiator2mm_D50cm_target6-15mm_gp/"
#rad_th = [2.0] #mm
#Dcm = [50] #cm
#conv_th = list(np.linspace(6.0, 15.0, 10)) #mm

#path = "/home/paterno/geant4-apps/PositronSource-build/output/results_6GeV_sigma0.5mm_W_radiator_var_D50cm_target_var_gp/"
#rad_th = [1.4, 1.6, 1.8, 2.0] #mm
#Dcm = [50] #cm
#conv_th = [9.0, 11.0, 13.0] #mm


useFullname = True


## Settings to read the built-in a BoxMesh scorer with Edep in the Absorber (Converter or Radiator)
if "crystalline" in path or "conventional" in path:
    tranvsizeX = 20      #volume tranversal X size (mm)
    tranvsizeY = 20      #volume tranversal Y size (mm)
else:
    tranvsizeX = 100     #volume tranversal X size (mm)
    tranvsizeY = 100     #volume tranversal Y size (mm)
#tranvsizeZ = conv_th    #volume tranversal Z size (mm) -> automatically retrieved below for each file!
IWantPlotVoxel = False
IWantPlotX0fract = False
X0 = 3.504               #crystal radiation length (mm)
normEvents = True        #normalize w.r.t. simulated events (set always to true to compare with previous results)
axis_proj = 2            #0: axial, 1:sagittal, 2:coronal
proj_type = 1            #1: mean, 2: sum
p2m = 1                  #voxels around which to mean
plot_proj = True
reader_verb = True
titstr = ''
cmap2D = 'viridis'
Wfig = 12
Hfig = 8
fs = 20
ms = 6

## Options for the plot of the positron distribution at the Converter exit
if "conventional" in path:
    inputType = "e-"
else:
    inputType = "all"    
bCutOnMomentum = False
MomentumCut = 5 #MeV
fitrange_ang = 17.45 #mrad (for angular distribution)
#if "crystalline" in path or "conventional" in path:
#    fitrange_ang = 10. #mrad (for angular distribution)    
ang_plot_lim = (-50, 150) #mrad
fitrange = 0.2 #cm (for spatial distribution)
nBins = 200 #for spatial distribution
xRange = (-0.45, 0.45) #cm
xRange1 = (-10.5, 10.5) #cm
yUpperLim = 14 #for profiles

## Settings for plotting Edep distribution
NbinsEdep = 50

## Settings fot plotting radius-vs-depth energy density
rLim = 1.5 #mm
bScaleWidth = False
bLog = False
upperLim = 20 #MeV/(mm^3 e-)
ind_r_lsit = [0, 7, 15, 22]

# Option to save the figures
saveFigs = True

# Option to convert file for RF-Track
convertFileToRFTrack = True
setGaussianTimeShape = True
addID = False
###############################################################################


# Define the output folder
outpath = path + "analysis_output/"
if not os.path.exists(outpath):
    os.makedirs(outpath)

# Quantity to be stored as results
case_list = []
yield_pos_list = []
yield_ele_list = []
yield_ph_list = []
yield_n_list = []
pos_mean_E_list = []
pos_spread_E_list = []
pos_mean_div_fit_list = []
pos_mean_size_fit_list = []
Edep_rad_list = []
Edep_conv_list = []
PEDD_list = []


# main loop
for itemi in rad_th:
    for itemj in Dcm:
        for itemk in conv_th:
            
            ## Read the root file with the scoring ntuples for virtual detectors
            if "crystalline" in path:
                if "mis" in path:
                    itemi_split = itemi.split('_')
                    name = "output_%sGeV_W%smm_crystalline_mis%srad_pot%sA" % (Ein, itemi_split[0], itemi_split[1], itemi_split[2])
                else:
                    if useFullname:
                        name = "output_%sGeV_W%smm_crystalline_mis0.000rad_pot0.050A" % (Ein, itemi)
                    else:
                        name = "output_%sGeV_W%smm_crystalline" % (Ein, itemi)         
            elif "conventional" in path:
                #name = "output_%sGeV_conventional_%s" % (Ein, itemi) #<--- used for the statistical analysis at optimal thickness
                name = "output_%sGeV_W%smm_conventional" % (Ein, itemi)
            else:
                name = "output_%sGeV_W%smm_D%scm_target%smm" % (Ein, round(itemi,1), round(itemj,0), round(itemk,1))
            purename = name.replace("output_", "")
            print("\nopening rf:", name ,".root file ...")

            rf = uproot.open(path + name + ".root")
            rf_content = [item.split(';')[0] for item in rf.keys()]
            df = rf['scoring_ntuple'].arrays(library='pd')

            if 'eventID' in df:
                Nevents = len(np.array(df.eventID.unique()))
            else:
                Nevents = 1e4
            print("Nevents:", Nevents) #number of simulated events (primary particles)

            df_rad_out = df[(df.screenID == 2)].copy().drop(["screenID"], axis=1)
            df_conv_in = df[(df.screenID == 0)].copy().drop(["screenID"], axis=1)
            df_conv_out = df[(df.screenID == 1)].copy().drop(["screenID"], axis=1)
            del df


            ##Plot the positron distribution at the Converter/Target exit 
            ##[and correlate it with the photon distribution at the exit of the radiator crystal]      
            
            ### Define data_in and data_out dataframes and add more parameters to them
            ### (adapt dataframes to use Mattia's scripts).        
            if df_rad_out.shape[0] > 0:
                data_in = df_rad_out.copy()
            else:
                if df_conv_out.shape[0] > 0:
                    data_in = df_conv_in.copy()
                else:
                    data_in = df_rad_out.copy() #set an empty dataframe, because this is a single-volume (crystalline or conventional) source  

            data_in["P"] = (data_in.px*data_in.px + data_in.py*data_in.py + data_in.pz*data_in.pz)**0.5 #MeV
            data_in = data_in[data_in.pz >= 0] #selecting only events (that should be) from the input beam
            if bCutOnMomentum:
                data_in = data_in[data_in.P >= MomentumCut] #selecting only events with momentum > MomentumCut
            data_in.pz = data_in.pz/1000 #MeV -> GeV
            data_in.x = data_in.x/10 #mm -> cm
            data_in.y = data_in.y/10 #mm -> cm
            data_in["pt"] = (data_in.px**2 + data_in.py**2)**0.5 #MeV
            data_in["thx"] = np.arctan(data_in.px/1000 / data_in.pz)*1000 #mrad
            data_in["thy"] = np.arctan(data_in.py/1000 / data_in.pz)*1000 #mrad
            data_in["tht"] = np.arctan(data_in.pt/1000 / data_in.pz)*1000 #mrad

            if df_conv_out.shape[0] > 0:
                data_out = df_conv_out.copy() 
            else: #this is the case of a single-volume (crystalline or conventional) source
                data_out = df_conv_in.copy()  

            data_out["P"] = (data_out.px*data_out.px+data_out.py*data_out.py+data_out.pz*data_out.pz)**0.5 #MeV
            if bCutOnMomentum:
                data_out = data_out[data_out.P >= MomentumCut] #selecting only events with momentum > MomentumCut
            data_out.pz = data_out.pz/1000 #MeV -> GeV
            data_out.x = data_out.x/10 #mm -> cm
            data_out.y = data_out.y/10 #mm -> cm
            data_out["pt"] = (data_out.px**2 + data_out.py**2)**0.5 #MeV
            data_out["thx"] = np.arctan(data_out.px/1000 / data_out.pz)*1000 #mrad
            data_out["thy"] = np.arctan(data_out.py/1000 / data_out.pz)*1000 #mrad
            data_out["tht"] = np.arctan(data_out.pt/1000 / data_out.pz)*1000 #mrad


            ### Plot the angular distribution            
            fig, ax = plt.subplots(nrows=2, sharex=True)

            divergenceIn = [0, 0]
            angleStdIn = [0, 0]
            divergence = [0, 0]
            angleStd = [0, 0]

            if len(data_in) > 0:
                series = [data_in[(data_in.particle=="gamma" if inputType=="all" else data_in.particle=="e-")].thx, 
                          data_in[(data_in.particle=="gamma" if inputType=="all" else data_in.particle=="e-")].thy]
                for i in (0, 1):
                    hist = ax[i].hist(series[i], bins=15000, range=(-200, 200), histtype="step", color="C0")
                    x, y = hist[1][:-1]+0.5*abs(hist[1][1]-hist[1][0]), hist[0]
                    x, y = x[(x>-fitrange_ang) & (x<fitrange_ang)], y[(x>-fitrange_ang) & (x<fitrange_ang)]
                    par, _ = curve_fit(sl.fGaus, x, y)
                    xplot = np.linspace(-fitrange_ang, fitrange_ang, 200)
                    xplot2 = np.linspace(-200, -fitrange_ang, 2000)
                    xplot3 = np.linspace(fitrange_ang, 200, 2000)
                    label = ("at crystal output"+r" ($\it{\gamma}$ only)" if inputType=="all" else \
                             "at amorphous target input"+r" ($\it{e}^-$ only)")+", Gaussian\n"+r"fit (in $\pm$ %.2f mrad)" % \
                             fitrange_ang+r" $\it{\sigma}$ = %.2f mrad" % abs(par[2])
                    ax[i].plot(xplot, sl.fGaus(xplot, *par), c="C0", alpha=0.5, label=label)
                    ax[i].plot(xplot2, sl.fGaus(xplot2, *par), c="C0", ls=":", alpha=0.5)
                    ax[i].plot(xplot3, sl.fGaus(xplot3, *par), c="C0", ls=":", alpha=0.5)
                    divergenceIn[i] = abs(par[2])
                    angleStdIn[i] = series[i].std()

                ax[i].set_yscale("linear")
                ax[i].set_xlim(ang_plot_lim)
                ax[i].set_title("                    "+("horizontal" if i==0 else "vertical"), loc="left")
                ax[i].legend(loc="upper right", fontsize="xx-small")
                ax[i].grid(True)

            series = [data_out[(data_out.particle=="e+")].thx, data_out[(data_out.particle=="e+")].thy]
            for i in (0, 1):
                hist = ax[i].hist(series[i], bins=800, range=(-1000, 1000), histtype="step", color="C1")
                x, y = hist[1][:-1]+0.5*abs(hist[1][1]-hist[1][0]), hist[0]
                x, y = x[(x>-fitrange_ang) & (x<fitrange_ang)], y[(x>-fitrange_ang) & (x<fitrange_ang)]
                par, _ = curve_fit(sl.fGaus, x, y)
                xplot = np.linspace(-fitrange_ang, fitrange_ang, 200)
                xplot2 = np.linspace(-200, -fitrange_ang, 2000)
                xplot3 = np.linspace(fitrange_ang, 200, 2000)
                label = "at amorphous target output"+r" ($\it{e}^+$ only)"+", Gaussian\n"+r"fit (in $\pm$ %.2f mrad)" % \
                        fitrange_ang+r" $\it{\sigma}$ = %.2f mrad" % abs(par[2])
                ax[i].plot(xplot, sl.fGaus(xplot, *par), c="C1", alpha=0.5, label=label)
                ax[i].plot(xplot2, sl.fGaus(xplot2, *par), c="C1", ls=":", alpha=0.5)
                ax[i].plot(xplot3, sl.fGaus(xplot3, *par), c="C1", ls=":", alpha=0.5)
                divergence[i] = abs(par[2])
                angleStd[i] = series[i].std()

                ax[i].set_yscale("linear")
                ax[i].set_xlim(ang_plot_lim)
                ax[i].set_title("                    "+("horizontal" if i==0 else "vertical"), loc="left")
                ax[i].legend(loc="upper right", fontsize="xx-small")
                ax[i].grid(True)

            fig.supylabel("Counts", fontsize="small")
            ax[1].set_xlabel(r"$\it{\theta}$ [mrad]")
            fig.tight_layout()

            if saveFigs:
                plt.savefig(outpath + 'eplusBeam_theta_distrib_' + purename + '.jpg')
            plt.close()

            print("e+ beam divergence along X and Y from Gauss fit in (-%.2f, %.2f) mrad range: [-%.2f, %.2f] mrad" \
                  % (fitrange_ang, fitrange_ang, divergence[0], divergence[1]))


            ### Plot the spatial distribution
            fig, ax = plt.subplots(ncols=2, nrows=3, figsize=(style['figure.figsize'][0], 14.0))
            fig.subplots_adjust(wspace=-1)

            beamSizeX, beamSizeY = [0, 0, 0], [0, 0, 0]

            ax[0, 0].set_title(r"$\it{e}^{+}$ only")
            ax[0, 1].set_title(r"$\it{e}^{+}$ only")

            ax[0, 1].hist2d(data_out[(data_out.particle=="e+")].x, data_out[(data_out.particle=="e+")].y, \
                            bins=nBins, range=(xRange, xRange), norm=LogNorm())
            ax[0, 1].grid(True)
            ax[0, 0].hist2d(data_out[(data_out.particle=="e+")].x, data_out[(data_out.particle=="e+")].y, \
                            bins=nBins, range=(xRange1, xRange1), norm=LogNorm())
            ax[0, 0].grid(True)

            gs = ax[1, 0].get_gridspec()
            for ax0 in ax[1, :]:
                ax0.remove()
            axdown = fig.add_subplot(gs[1, :])

            if len(data_in) > 0:
                hist = axdown.hist(data_in[(data_in.particle=="gamma" if inputType=="all" else data_in.particle=="e-")].x, \
                                   bins=nBins, histtype="step", range=xRange, density=True, color="C0")
                x, y = hist[1][:-1]+0.5*abs(hist[1][1]-hist[1][0]), hist[0]
                x, y = x[(x>-fitrange) & (x<fitrange)], y[(x>-fitrange) & (x<fitrange)]
                par, _ = curve_fit(sl.fGaus, x, y)
                xplot = np.linspace(-fitrange, fitrange, 200)
                xplot2 = np.linspace(xRange[0], -fitrange, 2000)
                xplot3 = np.linspace(fitrange, xRange[1], 2000)
                label = ("at crystal output"+r" ($\it{\gamma}$ only)" if inputType=="all" else \
                         "at amorphous target input"+r" ($\it{e}^-$ only)")+",\nGaussian fit "+r"$\it{\sigma}$ = %.3f mm" % abs(par[2]*10)
                axdown.plot(xplot, sl.fGaus(xplot, *par), c="C0", alpha=0.5, label=label)
                axdown.plot(xplot2, sl.fGaus(xplot2, *par), c="C0", ls=":", alpha=0.5)
                axdown.plot(xplot3, sl.fGaus(xplot3, *par), c="C0", ls=":", alpha=0.5)
                beamSizeX[0] = par[2]

            hist = axdown.hist(data_out[(data_out.particle=="e+")].x, bins=nBins, histtype="step", range=xRange, density=True, color="C1")
            x, y = hist[1][:-1]+0.5*abs(hist[1][1]-hist[1][0]), hist[0]
            x, y = x[(x>-fitrange) & (x<fitrange)], y[(x>-fitrange) & (x<fitrange)]
            par, _ = curve_fit(sl.fGaus, x, y)
            xplot = np.linspace(-fitrange, fitrange, 200)
            xplot2 = np.linspace(xRange[0], -fitrange, 2000)
            xplot3 = np.linspace(fitrange, xRange[1], 2000)
            label = "at amorphous target output"+r" ($\it{e}^{+}$ only)"+",\nGaussian fit "+r"$\it{\sigma}$ = %.3f mm" % abs(par[2]*10)
            axdown.plot(xplot, sl.fGaus(xplot, *par), c="C1", alpha=0.5, label=label)
            axdown.plot(xplot2, sl.fGaus(xplot2, *par), c="C1", ls=":", alpha=0.5)
            axdown.plot(xplot3, sl.fGaus(xplot3, *par), c="C1", ls=":", alpha=0.5)
            beamSizeX[1] = par[2]

            gs = ax[2, 0].get_gridspec()
            for ax0 in ax[2, :]:
                ax0.remove()
            axdown2 = fig.add_subplot(gs[2, :])

            if len(data_in) > 0:
                hist = axdown2.hist(data_in[(data_in.particle=="gamma" if inputType=="all" else data_in.particle=="e-")].y, \
                                    bins=nBins, histtype="step", range=xRange, density=True, color="C0")
                x, y = hist[1][:-1]+0.5*abs(hist[1][1]-hist[1][0]), hist[0]
                x, y = x[(x>-fitrange) & (x<fitrange)], y[(x>-fitrange) & (x<fitrange)]
                par, _ = curve_fit(sl.fGaus, x, y)
                xplot = np.linspace(-fitrange, fitrange, 200)
                xplot2 = np.linspace(xRange[0], -fitrange, 2000)
                xplot3 = np.linspace(fitrange, xRange[1], 2000)
                label = ("at crystal output"+r" ($\it{\gamma}$ only)" if inputType=="all" else \
                         "at amorphous target input"+r" ($\it{e}^-$ only)")+",\nGaussian fit "+r"$\it{\sigma}$ = %.3f mm" % abs(par[2]*10)
                axdown2.plot(xplot, sl.fGaus(xplot, *par), c="C0", alpha=0.5, label=label)
                axdown2.plot(xplot2, sl.fGaus(xplot2, *par), c="C0", ls=":", alpha=0.5)
                axdown2.plot(xplot3, sl.fGaus(xplot3, *par), c="C0", ls=":", alpha=0.5)
                beamSizeY[0] = par[2]

            hist = axdown2.hist(data_out[(data_out.particle=="e+")].y, bins=nBins, histtype="step", range=xRange, density=True, color="C1")
            x, y = hist[1][:-1]+0.5*abs(hist[1][1]-hist[1][0]), hist[0]
            x, y = x[(x>-fitrange) & (x<fitrange)], y[(x>-fitrange) & (x<fitrange)]
            par, _ = curve_fit(sl.fGaus, x, y)
            xplot = np.linspace(-fitrange, fitrange, 200)
            xplot2 = np.linspace(xRange[0], -fitrange, 2000)
            xplot3 = np.linspace(fitrange, xRange[1], 2000)
            label = "at amorphous target output"+r" ($\it{e}^{+}$ only)"+",\nGaussian fit "+r"$\it{\sigma}$ = %.3f mm" % abs(par[2]*10)
            axdown2.plot(xplot, sl.fGaus(xplot, *par), c="C1", alpha=0.5, label=label)
            axdown2.plot(xplot2, sl.fGaus(xplot2, *par), c="C1", ls=":", alpha=0.5)
            axdown2.plot(xplot3, sl.fGaus(xplot3, *par), c="C1", ls=":", alpha=0.5)
            beamSizeY[1] = par[2]

            ax[0, 0].set_ylabel(r"$\it{y}$ at amorphous"+"\ntarget output [cm]")
            ax[0, 0].set_xlabel(r"                 $\it{x}$ at amorphous target output [cm]", loc="left")
            axdown.set_xlabel(r"$\it{x}$ [cm]")
            axdown.set_ylabel(r"$\frac{\mathrm{d}\it{N}}{\it{N}\mathrm{d}\it{x}} \left[\frac{1}{\mathrm{cm}}\right]$")
            axdown.set_xlim((xRange[0], xRange[1]))
            axdown.set_ylim((0, yUpperLim))
            axdown.legend(loc="upper left", fontsize=16)
            axdown.grid(True)
            axdown2.set_xlabel(r"$\it{y}$ [cm]")
            axdown2.set_ylabel(r"$\frac{\mathrm{d}\it{N}}{\it{N}\mathrm{d}\it{y}} \left[\frac{1}{\mathrm{cm}}\right]$")
            axdown2.set_xlim((xRange[0], xRange[1]))
            axdown2.set_ylim((0, yUpperLim))
            axdown2.legend(loc="upper left", fontsize=16)
            axdown2.grid(True)
            fig.tight_layout()
            
            if saveFigs:
                plt.savefig(outpath + 'eplusBeam_spatial_distrib_' + purename + '.jpg')
            plt.close()

            beamSizeXStd = [data_in.x.std(), data_out[(data_out.particle=="e+")].x.std(), data_out[(data_out.particle=="e+")].x.std()]
            beamSizeYStd = [data_in.y.std(), data_out[(data_out.particle=="e+")].y.std(), data_out[(data_out.particle=="e+")].y.std()]

            print("e+ beam size along X and Y from Gauss fit in (-%.2f, %.2f) mm range: [-%.2f, %.2f] mm" \
                  % (fitrange, fitrange, np.abs(beamSizeX[1]*10), np.abs(beamSizeY[1]*10)))


            ## Plot the photon spectrum emitted by the Radiator crystal
            if len(data_in) > 0:
                df_ph_rad = data_in[(data_in.particle=="gamma")]
                nbin_gamma = 100
                range_gamma = (50, 2450)
                spectrum_Eph_rad, edges_Eph_rad = np.histogram(df_ph_rad.P, density=True, bins=nbin_gamma, range=range_gamma)
                bin_Eph_rad = edges_Eph_rad[:-1] + (edges_Eph_rad[1]-edges_Eph_rad[0])*0.5
                spectral_int_Eph_rad = spectrum_Eph_rad * bin_Eph_rad
                plt.figure(figsize=(8, 6))
                plt.plot(bin_Eph_rad, spectrum_Eph_rad, alpha=0.75, label='')
                #plt.legend()
                plt.xlabel('Energy [MeV]')
                plt.ylabel('1/N$\\times$dN/dE')
                #plt.ylabel('spectral intensity')
                plt.title('spectrum of the photons emitted by the radiator crystal')
                plt.grid(True)
                plt.yscale('log')
                if saveFigs:
                    plt.savefig(outpath + 'Eph_rad_' + purename + '.jpg')
                plt.close()
            
            
            ## Histograms of Edep per event in the Radiator and the Converter
            ## note: pay attention in case of read_PS...these plots make not much sense!

            ### Radiator
            if 'edep_rad' in rf_content:   
                df_edep_rad = rf['edep_rad'].arrays(library='pd')
                edep_rad = np.array(df_edep_rad["edep"])
                print("total Edep in the Radiator per event:", round(np.sum(edep_rad)/Nevents, 2), "MeV/e-")
                if len(df_edep_rad) > 0:
                    edep_rad_cut = edep_rad[edep_rad > 1]
                    plt.figure(figsize=(8, 6))
                    h_rad = plt.hist(edep_rad_cut, bins=NbinsEdep,
                                     histtype="step", label="", density=False,
                                     edgecolor="blue", linewidth=1.2, alpha=1.)
                    plt.title('Edep per event in the radiator', fontsize=fs)
                    plt.xlabel("Edep [MeV]", fontsize=fs)
                    plt.ylabel("Counts", fontsize=fs)
                    plt.xticks(fontsize=fs, rotation=0)
                    plt.yticks(fontsize=fs, rotation=0)
                    plt.grid(which="major", color="gray", linestyle="--", linewidth=1)
                    if saveFigs:
                        plt.savefig(outpath + 'Edep_distrib_rad_' + purename + '.jpg')
                    plt.close()

            ### Converter
            if 'edep_conv' in rf_content:
                df_edep_conv = rf['edep_conv'].arrays(library='pd')
                edep_conv = np.array(df_edep_conv["edep"])
                print("total Edep in the Converter per event:", round(np.sum(edep_conv)*1.e-3/Nevents, 2), "GeV/e-")
                if len(df_edep_conv) > 0 and not "_GT" in name:       
                    edep_conv_cut = edep_conv[edep_conv > 0]
                    plt.figure(figsize=(8, 6))
                    h_conv = plt.hist(edep_conv_cut, bins=NbinsEdep,
                                      histtype="step", label="", density=False,
                                      edgecolor="blue", linewidth=1.2, alpha=1.)
                    plt.title('Edep per event in the converter', fontsize=fs)
                    plt.xlabel("Edep [MeV]", fontsize=fs)
                    plt.ylabel("Counts", fontsize=fs)
                    plt.xticks(fontsize=fs, rotation=0)
                    plt.yticks(fontsize=fs, rotation=0)
                    plt.grid(which="major", color="gray", linestyle="--", linewidth=1)
                    if saveFigs:
                        plt.savefig(outpath + 'Edep_distrib_conv_' + purename + '.jpg')
                    plt.close()
                    

            ### Read the built-in BoxMesh scorer with Edep in the Absorber
            
            # read BoxMesh scorer
            if itemk != '':                
                tranvsizeZ = itemk
            else:
                if type(itemi) == str:
                    tranvsizeZ = np.double(itemi.split('_')[0])
                else:   
                    tranvsizeZ = itemi
            filename = path + name.replace("output", "Edep") + ".txt"
            results_BoxMesh = read_Edep_BoxMesh(filename, normEvents, Nevents, 
                                                tranvsizeX, tranvsizeY, tranvsizeZ)
            df_BoxMesh = results_BoxMesh[0]
            voxel_coord = results_BoxMesh[1]
            del results_BoxMesh

            # retrieve variables form BoxMesh scorer
            eDep = df_BoxMesh.eDep
            eDepDensity = df_BoxMesh.eDepDensity
            x = voxel_coord[0]
            y = voxel_coord[1]
            z = voxel_coord[2]
            if x.shape[0] > 1:
                dx = x[1] - x[0]
            else:
                dx = 0
            if y.shape[0] > 1:
                dy = y[1] - y[0]
            else:
                dy = 0
            if z.shape[0] > 1:
                dz = z[1] - z[0]
            else:
                dz = 0

            # variables for contour plot of Edep
            nx, ny , nz = [voxel_coord[i].shape[0] for i in range(3)]
            X, Z = np.meshgrid(x, z)
            icx = 0 if nx == 1 else int(X.shape[1]/2+1)
            eDepMap = np.transpose(np.array(eDep).reshape(nx,ny,nz)[icx,:,:])
            eDepDensityMap = np.transpose(np.array(eDepDensity).reshape(nx,ny,nz)[icx,:,:])

            # Zprofile of Edep
            ic = int(Z.shape[1]/2+1)
            profileZ_eDep = np.zeros(nz)
            profileZ_eDepDensity = np.zeros(nz)
            for j in range(nz):
                if proj_type == 1:
                    profileZ_eDep[j] = np.mean(eDepMap[j][ic-p2m:ic+p2m])
                    profileZ_eDepDensity[j] = np.mean(eDepDensityMap[j][ic-p2m:ic+p2m])
                else:
                    profileZ_eDep[j] = np.sum(eDepMap[j][ic-p2m:ic+p2m])
                    profileZ_eDepDensity[j] = np.sum(eDepDensityMap[j][ic-p2m:ic+p2m])

            PEDD = np.max(np.array(eDepDensity))
            print("PEDD:", round(PEDD,2), "MeV/(mm^3*e-)")


            ### Contour and plofile plot of Edep (density) distribution
            
            # set options
            if normEvents:
                clbl = "Edep density per event [MeV/(mm$^3$ e-]"
            else:
                clbl = "Edep density [MeV/mm$^3$]"

            # contour plot of Edep
            if eDepDensityMap.shape[0]>1 and eDepDensityMap.shape[1]>1: 
                plt.figure(figsize=(Wfig, Hfig))
                Nlev = 50
                if IWantPlotVoxel:
                    plt.contourf(eDepDensityMap, Nlev, cmap='viridis')
                    xlabel = "X [voxel]"
                    ylabel = "Z [voxel]"        
                else:
                    plt.contourf(X, Z, eDepDensityMap, Nlev, cmap='viridis')
                    xlabel = "X [mm]"
                    ylabel = "Z [mm]"
                #plt.gca().set_aspect('equal')
                plt.gca().invert_yaxis()
                cbar = plt.colorbar()
                cbar.ax.tick_params(labelsize=fs)
                cbar.set_label(clbl, fontsize=fs, rotation=90)
                plt.title(titstr, fontsize=fs)
                plt.xlabel(xlabel, fontsize=fs)
                plt.ylabel(ylabel, fontsize=fs, wrap=True)
                plt.xticks(fontsize=fs, rotation=0)
                plt.yticks(fontsize=fs, rotation=0)
                plt.gca().xaxis.set_minor_locator(AutoMinorLocator(5))
                plt.gca().yaxis.set_minor_locator(AutoMinorLocator(5))
                plt.gca().tick_params(axis="both", which='major', direction='in', length=8)
                plt.gca().tick_params(axis="both", which='minor', direction='in', length=4)
                plt.grid(False)
                if saveFigs:
                    plt.savefig(outpath + 'Edep_proj_' + purename + '.jpg')
                plt.close()

            # plot Zprofile of Edep
            if IWantPlotVoxel:
                xplot = np.linspace(0, nz-1, nz)
                ylabel = "Z [voxel]"  
            else:
                xplot = z
                ylabel = "Z [mm]"
            plt.figure(figsize=(Wfig, Hfig))
            plt.plot(xplot, profileZ_eDepDensity, label="",
                     color="blue", linestyle='-', linewidth=2, alpha=0.75,
                     marker='', markersize=ms, markerfacecolor="Blue")
            plt.title(titstr, fontsize=fs)
            plt.xlabel(ylabel, fontsize=fs, wrap=True)
            plt.ylabel(clbl, fontsize=fs, wrap=True)
            plt.xticks(fontsize=fs, rotation=0)
            plt.yticks(fontsize=fs, rotation=0)
            plt.gca().xaxis.set_minor_locator(AutoMinorLocator(5))
            plt.gca().yaxis.set_minor_locator(AutoMinorLocator(5))
            plt.gca().tick_params(axis="both", which='major', direction='in', length=8)
            plt.gca().tick_params(axis="both", which='minor', direction='in', length=4)
            plt.grid(True)
            if saveFigs:
                plt.savefig(outpath + 'Edep_profileZ_' + purename + '.jpg')
            plt.close()


            ### radius-vs-depth energy density

            # radius-vs-depth energy density
            bData = df_BoxMesh["r"] < rLim if rLim != None else df_BoxMesh["r"] < 1e5
            Z, R = np.meshgrid(
                sorted(list(set(df_BoxMesh[bData]["z"]))),
                sorted(list(set(df_BoxMesh[bData]["r"]))),
            )
            ED = interpolate.griddata(
                (df_BoxMesh[bData]["z"], df_BoxMesh[bData]["r"]),
                 df_BoxMesh[bData]["eDepDensity"], (Z, R), method='nearest'
            )
            
            if ED.shape[0]>1 and ED.shape[1]>1: 
                width = 5*max(df_BoxMesh["z"]) if bScaleWidth else 10
                fig, ax = plt.subplots(num="mesh", figsize=[width, 8], nrows=2, sharex=True, constrained_layout=True)
    
                # 2d plot
                plotmesh = ax[0].pcolormesh(Z, R, ED, shading="nearest", norm=LogNorm() if bLog else None)
                cbar = plt.colorbar(plotmesh, ax=ax[0])
                ax[0].grid()
                ax[0].set_ylim((0, max(df_BoxMesh[bData]["r"])))

                # 1d projection - energy density vs depth
                dr = (dx**2 + dy**2)**0.5
                for ind_r in ind_r_lsit:
                    ind_r_eff = int(ind_r * 0.1 / dr)
                    r_eff = dr * ind_r_eff
                    bData = (df_BoxMesh["r"] >= (dr * (ind_r_eff - 0.5))) & \
                            (df_BoxMesh["r"] <= (dr * (ind_r_eff + 0.5)))   
                    temp = df_BoxMesh[bData][["z", "eDepDensity"]].sort_values(by=["z"])
                    zlist = np.array(temp.z.unique())
                    ymean = [np.mean(np.array(temp[temp["z"] == zi]["eDepDensity"])) for zi in zlist]
                    ax[1].plot(
                        #df_BoxMesh[bData]["z"], df_BoxMesh[bData]["eDepDensity"],
                        zlist, ymean, #corrected by gpaterno
                        label = "r = %.2f mm" % r_eff,
                        marker='', linestyle='-'
                    )
                #ax[1].set_ylim((0, upperLim))
                ax[1].grid()
                ax[1].legend(loc="upper left", fontsize=14)

                fs_local = 18
                ax[1].set_xlabel("z [mm]", fontsize=fs_local)
                ax[0].set_ylabel("r [mm]", fontsize=fs_local)
                ax[1].set_ylabel(r"eDep [MeV/(mm$^3$ e-)]", fontsize=fs_local)
                cbar.set_label(r"eDep [MeV/(mm$^3$ e-)]", rotation=90, labelpad=16, fontsize=fs_local)
                #fig.suptitle("e- beam features: %.2f GeV, %.1f mm" % (Ein, sigma_in))
                if saveFigs:
                    plt.savefig(outpath + 'Edep_profileR_' + purename + '.jpg')
                plt.close()


            ## Store the positron beam phase-space in RF-Tack format and prepare the variables to save
            
            # Get the positrons
            positrons = data_out[data_out.particle=="e+"].copy()
            print("yield_e+:", round(positrons.shape[0] / Nevents, 2))
           
        
            # Calculate positron beam energy distribution features
            E_pos = (positrons.P**2 + 0.511**2)**0.5
            Emean_pos = np.mean(E_pos)
            Estdev = np.std(E_pos)
            Espread_pos = Estdev / Emean_pos
            print("Emean_pos:", round(Emean_pos, 2), "MeV")
            print("Espread_pos:", round(Espread_pos, 2), "(sigma/mu)")
            angle_cut = np.mean(divergence) #mrad   
            #angle_cut = 140 #mrad
            positrons_cut = positrons[positrons["tht"] < angle_cut]
            E_pos_cut = (positrons_cut.P**2 + 0.511**2)**0.5 #MeV
            Emean_pos_cut = np.mean(E_pos_cut)
            Estdev_cut = np.std(E_pos_cut)
            Espread_pos_cut = Estdev_cut / Emean_pos_cut
            print("Emean_pos_cut (tht < %.2f mrad):" % angle_cut, round(Emean_pos_cut, 2), "MeV")
            print("Espread_pos_cut (tht < %.2f mrad):" % angle_cut, round(Espread_pos_cut, 2), "(sigma/mu)")
            Eth = 100 # MeV
            NposEth = len(E_pos.values[E_pos.values < Eth]) / len(E_pos.values)
            NsigmaEth = Eth / (Emean_pos + Estdev)
            print("fraction of positrons with energy < %.0f MeV (%.2f*sigma): %.2f\n" % (Eth, NsigmaEth, NposEth))

            # Plot histogram of positrons energy
            nbin_pos = 100
            range_pos = (0, 100)
            IWantDensity = False
            plt.figure(figsize=(9, 6))
            h = plt.hist(E_pos, density=IWantDensity, bins=nbin_pos, range=range_pos, alpha=0.5, \
                         label='positron spectrum')
            h2 = plt.hist(E_pos_cut, density=IWantDensity, bins=nbin_pos, range=range_pos, alpha=0.5, \
                          label='positron spectrum within %.2f mrad' % (angle_cut))
            plt.legend()
            plt.xlabel('Energy [MeV]')
            plt.ylabel('Counts [arb. units]')
            plt.title('')
            plt.grid(True)
            plt.yscale('log')
            if saveFigs:
                plt.savefig(outpath + 'eplusBeam_E_distrib_' + purename + '.jpg')
            plt.show()

            if not IWantDensity:
                print("sum of positron spectrum: %d" % sum(h[0]))
                print("sum of positron spectrum within %.2f mrad: %d\n" % (angle_cut, sum(h2[0])))
            else:
                print("\n")
            """
            Eedges, spectrum = h[1], h[0]
            output_file = outpath + 'positron_spectrum_' + purename + '.txt'
            write_spectrum_to_G4file(Eedges, spectrum, output_file)
            _, spectrum_cut = h2[1], h2[0]
            output_file_cut = outpath + 'positron_spectrum_cut_%.1fmrad_' % angle_cut + purename + '.txt'
            write_spectrum_to_G4file(Eedges, spectrum_cut, output_file_cut)
            """
            
            # Prepare the file to feed the RF-Track code to tack the positrons inside the Capture System/Positron Linac with
            """
            NOTE: In RF-Track x'=px/pz=tan(thx) and y'=py/pz=tan(thy) (see the manual at page 22).
                  The approximation x'=thx and y'=thy holds for small-angle only,
                  which is the case of particle beams in accelators. However, 
                  for example in a drift, we correctly have x2=x1+x'*(z2-z1) and y2=y1+y'*(z2-z1).
                  In our case, positrons can exit the crystal also at large angles 
                  and, for some of them, px and/or py can be (much) greater than pz. Thus,
                  x' and y' are not angles, but tangents (which can tend to infinity) and 
                  the fact that they are expressed in mrad means that they are multiplied by 1e-3.
                  These particles would be lost in the capture section. However, if we want to
                  calculate the tranverse phase space of the positron beam at the exit of the crystal,
                  we must calculate properly (with arctan) thx, thy and set a proper angular cut.
            """
            if convertFileToRFTrack:    
                ## Add variables to the original positron dataframe
                positrons["xp[mrad]"] = (positrons.px/(positrons.pz*1000))*1000
                positrons["yp[mrad]"] = (positrons.py/(positrons.pz*1000))*1000
                positrons["p[MeV/c]"] = positrons.P
                positrons["#x[mm]"] = positrons["x"]*10 #this is because they were previously converted form [mm] to [cm]
                positrons["y[mm]"] = positrons["y"]*10

                ## Select only a set of variables
                if addID:
                    positrons["ID"] = list(range(1, len(positrons)+1))
                    positrons_out = positrons[["#x[mm]", "xp[mrad]", "y[mm]", "yp[mrad]", "p[MeV/c]", "ID"]]
                else:
                    positrons_out = positrons[["#x[mm]", "xp[mrad]", "y[mm]", "yp[mrad]", "p[MeV/c]"]]

                ## Guassian sampling in time, to keep the RF phases of the Hybrid similar to the conventional scheme
                c = 299792458 #m/s
                if setGaussianTimeShape:
                    mean_value = 17.762 #mm/c
                    std_value = 1.1950 #mm/c
                    num_sample = len(positrons_out)
                    gaussian_samples = np.random.normal(mean_value, std_value, num_sample)
                    t = gaussian_samples
                else:
                    t = positrons.t * 1e-6 * c #ns -> mm/c
                positrons_out.insert(4, 't[mm/c]', t)
                
                ## Save the positron phase-space to a txt file
                name_out = name.replace('output', 'positrons')
                positrons_out.to_csv(outpath + name_out + ".dat", index=False, sep=' ') 
                print("positron beam phase-space converted to RF-Track format!")
            

            # Plot the positron beam time distribution
            plt.figure(figsize=(10, 8))
            trLim = (-3*np.std(positrons.t) + np.mean(positrons.t), 3*np.std(positrons.t) + np.mean(positrons.t))
            plt.hist(positrons.t, density=False, bins=100, range=trLim, label='positrons original time distribution')
            if convertFileToRFTrack and setGaussianTimeShape:
                plt.hist(positrons_out["t[mm/c]"]*1e6/c, density=False, bins=100, label='positrons distribution after Gaussian shaping')
            plt.legend(fontsize=12)
            plt.xlabel('t [ns]')
            plt.ylabel('Counts [arb. units]')
            plt.title('')
            plt.grid(True)
            #plt.yscale('log')
            if saveFigs:
                plt.savefig(outpath + 'eplusBeam_t_distrib_' + purename + '.jpg')
            plt.close()
            
            
            # Fill the result list with those of the current simulation
            case_list.append(purename)
            yield_pos_list.append(round(data_out[(data_out.particle=="e+")].shape[0] / Nevents, 2))
            yield_ele_list.append(round(data_out[(data_out.particle=="e-")].shape[0] / Nevents, 2))
            yield_ph_list.append(round(data_out[(data_out.particle=="gamma")].shape[0] / Nevents, 2))
            yield_n_list.append(round(data_out[(data_out.particle=="neutron")].shape[0] / Nevents, 2))
            pos_mean_E_list.append(round(Emean_pos, 2))
            pos_spread_E_list.append(round(Espread_pos, 2))   
            pos_mean_div_fit_list.append(round(np.mean(divergence), 2))
            pos_mean_size_fit_list.append(round(np.mean([np.abs(beamSizeX[1]*10), np.abs(beamSizeY[1]*10)]), 2))
            Edep_rad_list.append(round(np.sum(edep_rad)/Nevents, 2))
            Edep_conv_list.append(round(np.sum(edep_conv)*1e-3/Nevents, 2))
            PEDD_list.append(round(PEDD, 2))

            
## Store the results to a dictionary
results = { "case": case_list,
            "yield_e+": yield_pos_list,
            "yield_e-": yield_ele_list,
            "yield_ph": yield_ph_list,
            "yield_n": yield_n_list,
            "e+_mean_E[MeV]": pos_mean_E_list,
            "e+_spread_E[sigma/mu]": pos_spread_E_list,
            "e+_mean_div_fit[mrad]": pos_mean_div_fit_list,
            "e+_mean_size_fit[mm]": pos_mean_size_fit_list,
            "Edep_rad[MeV/e-]": Edep_rad_list,
            "Edep_conv[GeV/e-]": Edep_conv_list,
            "PEDD[MeV/(mm^3*e-)]": PEDD_list
          }

# Get a dataframe of the results
df_results = pd.DataFrame.from_dict(results)
df_results

# plot some variables in the dataframe
if len(case_list[0]) > 10:
    if len(rad_th) > 1:
        case_list_short = [item.split('W')[1].split('_')[0] for item in case_list]
    elif len(Dcm) > 1:
        case_list_short = [item.split('_D')[1].split('_')[0] for item in case_list]
    elif len(conv_th) > 1:
        case_list_short = [item.split('target')[1] for item in case_list]
    else:
        case_list_short = range(len(case_list))
else:
    case_list_short = case_list
rot_value = 45
from matplotlib.ticker import FormatStrFormatter
fig, ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel('case')
ax1.set_ylabel('e+ yield', color=color)
ax1.plot(yield_pos_list, color=color, 
         linestyle='-', linewidth=2, marker='o', markersize=10)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_xticks(range(len(case_list_short)))
ax1.set_xticklabels(case_list_short, rotation=rot_value)
ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
color = 'tab:blue'
ax2 = ax1.twinx()
ax2.set_xticks(range(len(case_list_short)))
ax2.set_xticklabels(case_list_short, rotation=rot_value)
ax2.set_ylabel('PEDD [MeV/(mm$^3$ e-)]', color=color)
ax2.plot(PEDD_list, color=color, 
         linestyle='-', linewidth=2, marker='o', markersize=10)
ax2.tick_params(axis='y', labelcolor=color)
fig.tight_layout()
plt.title("e- beam features: %.2f GeV, %.1f mm" % (Ein, sigma_in))
if saveFigs:
    plt.savefig(outpath + 'yield_and_PEDD_vs_case' + '.jpg')
plt.close()
fig, ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel('case')
ax1.set_ylabel('e+ size [mm]', color=color)
ax1.plot(pos_mean_size_fit_list, color=color, 
         linestyle='-', linewidth=2, marker='o', markersize=10)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_xticks(range(len(case_list_short)))
ax1.set_xticklabels(case_list_short, rotation=rot_value)
ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
color = 'tab:blue'
ax2 = ax1.twinx()
ax2.set_xticks(range(len(case_list_short)))
ax2.set_xticklabels(case_list_short, rotation=rot_value)
ax2.set_ylabel('Edep [GeV/e-]', color=color)
if np.array(Edep_conv_list).any() == 0:
    var_to_plot = np.array(Edep_rad_list)*1e-3
else:
    var_to_plot = np.array(Edep_conv_list)
ax2.plot(var_to_plot, color=color, 
         linestyle='-', linewidth=2, marker='o', markersize=10)
ax2.tick_params(axis='y', labelcolor=color)
fig.tight_layout()
plt.title("e- beam features: %.2f GeV, %.1f mm" % (Ein, sigma_in))
if saveFigs:
    plt.savefig(outpath + 'eplusSize_and_EdepConv_vs_case' + '.jpg')
plt.close()
fig, ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel('case')
ax1.set_ylabel('e+ yield', color=color)
ax1.plot(yield_pos_list, color=color, 
         linestyle='-', linewidth=2, marker='o', markersize=10)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_xticks(range(len(case_list_short)))
ax1.set_xticklabels(case_list_short, rotation=rot_value)
ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
color = 'tab:blue'
ax2 = ax1.twinx()
ax2.set_xticks(range(len(case_list_short)))
ax2.set_xticklabels(case_list_short, rotation=rot_value)
ax2.set_ylabel('Edep [GeV/e-]', color=color)
if np.array(Edep_conv_list).any() == 0:
    var_to_plot = np.array(Edep_rad_list)*1e-3
else:
    var_to_plot = np.array(Edep_conv_list)
ax2.plot(var_to_plot, color=color, 
         linestyle='-', linewidth=2, marker='o', markersize=10)
ax2.tick_params(axis='y', labelcolor=color)
fig.tight_layout()
plt.title("e- beam features: %.2f GeV, %.1f mm" % (Ein, sigma_in))
if saveFigs:
    plt.savefig(outpath + 'yield_and_Edep_vs_case' + '.jpg')
plt.close()


# Export the results in a json file.
# pickle works pretty much the same, 
# but it I cannot directly have a look at it.
import json 
print("\nWriting json file with these data:")
for s in results:
    print("\t", s, ":", results[s])
outname = outpath + "results"                  
with open(outname + ".json", "w") as file:    
    json.dump(results, file)
file.close()

"""
# Print the dataframe
import df2img
fig = df2img.plot_dataframe(
    df_results.round(2),
    row_fill_color=("#ffffff", "#d7d8d6"),
    fig_size=(1500, 500),
)
df2img.save_dataframe(fig=fig, filename=outpath+"df_results.jpg")
"""


# toc
print('\nThe execution took', round(time.time()-start, 0), 'seconds.\n')
