#######################################################################################################
####### Set of functions useful to analyse Geant4 simulations #########################################
####### Author: Gianfranco Paternò (paterno@fe.infn.it), last update: 23/12/2025 ######################
#######################################################################################################


############# Set of functions to elaborate numpy arrays and lists ####################################
def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation of an array.
    They weights are in effect first normalized so that they 
    sum to 1 (and so they must not all be 0).
    values, weights -- NumPy ndarrays with the same shape.
    The calculation is fast and numerically precise.
    """
    import numpy as np
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)
    return (average, np.sqrt(variance))
    

def print_array_size(arr):
    """
    Function to print the size [MB] of a list.
    """
    import numpy as np 
    x = np.array(arr)
    print("Size of the array = number of elements:", x.size)
    print("Memory size of one array element in bytes:", x.itemsize)
    print("Memory size of numpy array in MB:", round(x.size * x.itemsize / 1048576, 2))
    
    
def list_flatten(myList):
    """
    Function to flatten a list of lists
    """
    flat_list = [item for sublist in myList for item in sublist]
    return flat_list


def find(cond, N=1e9):
    """
    Function equivalent to Matlab's one that returns
    a list with the indices of an array/matrix
    that satisfy a given condition.
    """
    import numpy as np 
    cond = cond.reshape(np.size(cond), 1).ravel()
    index = list(np.argwhere(cond).ravel())
    if N < len(index):
        return index[:N]
    else:
        return index


def find_list_max(myList):
    """
    Function to find max and imax of a list.
    """ 
    the_max = myList[0]
    imax = 0
    for i in range(1, len(myList)):
        if myList[i] > the_max:
            the_max = myList[i]
            imax = i
    return the_max, imax


def Average(lst):
    """
    Function to calculate the average of a list.
    """ 
    return sum(lst) / len(lst)


def from_edges_to_bins(edges): 
    """
    Function to calculate bin centers from bin edges.
    """
    import numpy as np    
    # Sanity Check
    assert (
        type(edges) == np.ndarray and edges.ndim == 1 and edges.shape[0] > 1
    ), "the input parameter must be a 1D numpy array with size > 1"  
    # Calculation of bins
    bin_width = edges[1] - edges[0]
    bins = edges + bin_width*0.5
    bins = bins[:-1]
    return bins

 
def rebin(x0, y0, nBin):
    """
    Function to rebin an histogram (downsampling only),
    nBin is downsampling ratio.
    """
    import numpy as np
    import math    
    x, y = np.empty(0), np.empty(0)
    for i in range(0, math.floor(len(x0)/nBin)*nBin, nBin):
        x = np.append(x, sum(x0[i:(i+nBin)])/nBin)
        y = np.append(y, sum(y0[i:(i+nBin)]))
    return x, y


def histc(x, edges):
    """
    Function to get the number of counts within 
    two consecutive bin edges. It is similar to that of Matlab.
    pyplot.hist and numpy.histogram functions do the same.
    """
    import numpy as np
    map_to_bins = np.digitize(x, edges)
    # Get indices of the bins to which each value in input array belongs
    res = np.zeros(edges.shape)
    for el in map_to_bins:
        res[el-1] += 1 # Increment appropriate bin
    return res[:-1], map_to_bins


def array_shift(my_array, th):
    """
    Function to shift a numpy array
    based on a threshold value
    """
    import numpy as np
    N = my_array.shape[0]
    shifted_array = np.zeros((N,))
    jth = 0
    for j in range(N):
        if my_array[j] > th:
            jth = j
            break
    for j in range(N):
        shifted_array[j-jth] = my_array[j]
    return shifted_array
    

def TProfile2D(X, Y, NbinX, Xmin, Xmax, NbinY, Ymin, Ymax, Z):
    """
    My algorithm to mimic the root TPofile2D object.
    It is similar to the one implemented in Matlab.
    """
    import numpy as np       
    if (len(X) != len(Y) or len(X) != len(Z)):
        print('X, Y and Z must have the same length')
    N = len(X)    
    def histc(x, edges):
        map_to_bins = np.digitize(x, edges)
        res = np.zeros(edges.shape)
        for el in map_to_bins:
            res[el-1] += 1
        return res[:-1], map_to_bins   
    XbinEdges = np.linspace(Xmin, Xmax, NbinX)
    YbinEdges = np.linspace(Ymin, Ymax, NbinY)  
    [nx, idx] = histc(X, XbinEdges)
    [ny, idy] = histc(Y, YbinEdges)   
    countMatrix = np.zeros([NbinY, NbinX]) 
    Zmatrix = np.zeros([NbinY, NbinX]) 
    for q in range(N):                                
        if (idx[q] != 0 and idy[q] != 0):
            i = idy[q] - 1
            j = idx[q] - 1 
            #print(q, i, j)
            countMatrix[i,j] = countMatrix[i,j] + 1
            Zmatrix[i,j] = Zmatrix[i,j] + Z[q]
    mask = np.squeeze(np.array([countMatrix > 0]))
    countMatrix = np.ma.MaskedArray(countMatrix, ~mask)
    Zmatrix_norm = Zmatrix / countMatrix
    return Zmatrix, Zmatrix_norm, XbinEdges, YbinEdges


def elliptical_selection(x, y, g_ell_center, g_ell_width, g_ell_height, angle, ax=0):
    """
    Custom function to select points with (x,y) coordinates within a given ellipse
    defined by its center, size of axes, and anti-clockwise tilt angle [deg].
    """
    import numpy as np
    import matplotlib.patches as patches
    cos_angle = np.cos(np.radians(angle))
    sin_angle = np.sin(np.radians(angle))
    xc = x - g_ell_center[0]
    yc = y - g_ell_center[1]
    xct = xc * cos_angle - yc * sin_angle
    yct = xc * sin_angle + yc * cos_angle 
    rad_cc = (xct**2/(g_ell_width/2.)**2) + (yct**2/(g_ell_height/2.)**2)
    mask = rad_cc <= 1
    xin = x[mask]
    yin = y[mask]
    g_ellipse = patches.Ellipse(g_ell_center, g_ell_width, \
                                g_ell_height, angle=angle, \
                                fill=False, edgecolor='red', linewidth=2)
    if ax:
        ax.add_patch(g_ellipse)
    return xin, yin, mask
#######################################################################################################

############# Set of functions present also in utils_testbeam #########################################
def smooth(y, box_pts):
    """
    Function to smooth data (y).
    box_pts acts a smooth factor.
    """
    import numpy as np
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def unique_with_global_threshold(arr, threshold):
    """
    It takes a list and returns the values
    that are unique within a given threshold.
    """
    import numpy as np
    arr_sorted = np.sort(arr)
    unique_values = []
    # Initialize a temporary group with the first value
    group = [arr_sorted[0]]
    # Loop
    for value in arr_sorted[1:]:
        # Compare the new value with the current group
        if abs(value - group[0]) <= threshold:
            group.append(value)
        else:
            # If the value is too far, append the mean of the group
            # to unique_values and start a new group
            unique_values.append(np.mean(group))
            group = [value]
    # Append the last group
    unique_values.append(np.mean(group)) 
    # Return
    return np.array(unique_values)


def projectDistZ3(x1, y1, z1, x2, y2, z2, z3):
    """
    General version of projectDistZ.
    Developed by gpaterno.
    """
    import numpy as np
    mx = (x2-x1)/(z2-z1)
    x3 = x1 + mx*z3    
    my = (y2-y1)/(z2-z1)
    y3 = y1 + my*z3
    return (x3, y3)


def get_theta_angles(x1, y1, x2, y2, d12):
    """
    Developed by gpaterno.
    """
    import numpy as np
    thetaX = np.arctan((x2-x1)/d12) #rad
    thetaY = np.arctan((y2-y1)/d12) #rad
    return (thetaX, thetaY)


def myGauss(x, a, mu, sigma):
    """Simple Gaussian function."""
    import numpy as np
    fun = a*np.exp(-(x-mu)**2/(2*sigma**2))
    return fun


def myLandau(x, a, b, c):
    """
    See: https://en.wikipedia.org/wiki/Landau_distribution.
    Note (gpaterno): if c<0 => left tail, if => c>0 right tail (True Laundau).
    """
    import numpy as np
    t = (x - b) / c
    fun = np.where(t < -5., 0., a*np.sqrt(1/(2*np.pi))*np.exp(-0.5*(t+np.exp(-t))))
    return fun


def myRayleigh(x, a, b):
    """
    Used in muonsPropagation to fit r=sqrt(x^2+y^2).
    """
    from scipy.stats import rayleigh
    fun = a*rayleigh.pdf(x, loc=0, scale=b)
    return fun
#######################################################################################################

############# Set of functions to do useful and most used plots #######################################
def simple_plot(x, y, yerr=None, icludeErrors=False, \
                xlbl='x', ylbl='y', fs=14, mytitle='', showPlot=True, \
                myFigSize=(8, 5), saveFig=False, figname='mySimplePlot'):
    """
    Function to plot y as a function of x. Errors can be also included.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator
    # Plot
    fig = plt.figure(figsize=myFigSize)
    if icludeErrors and yerr.size>0:
        plt.errorbar(x, y, yerr, \
                     linestyle='', linewidth=1, color='b', \
                     marker='o', ms=4, ecolor='b', capsize=0.0)
    else:
        plt.plot(x, y, \
                 color='b', linestyle='-', linewidth=2, \
                 marker='', markersize=8,  markerfacecolor='b')
    plt.title(mytitle, fontsize=fs)
    plt.xlabel(xlbl, fontsize=fs)
    plt.ylabel(ylbl, fontsize=fs)
    plt.xticks(fontsize=fs, rotation=0)
    plt.yticks(fontsize=fs, rotation=0)
    plt.gca().xaxis.set_minor_locator(AutoMinorLocator(5))
    plt.gca().yaxis.set_minor_locator(AutoMinorLocator(5))
    plt.gca().tick_params(axis="both", which='major', direction='in', length=8)
    plt.gca().tick_params(axis="both", which='minor', direction='in', length=4)
    plt.grid(which="major", color="gray", linestyle="--", linewidth=1)
    if saveFig:
        plt.savefig(figname+".jpg")
    ax = plt.gca()
    if showPlot:
        plt.show()
    #Return
    return ax
    

def simple_hist(values, Nbins=100, myrange=None, IWantDensity=False, \
                xlbl='x', fs=14, bw=0.5, mylabel='data', mytitle='', showPlot=True, \
                icludeErrors=False, IWantGaussianFit=False, plotFitPar=False, \
                myFigSize=(8, 5), saveFig=False, figname='mySimpleHist'):
    """
    Function to plot the histogram of passed values. A Gaussian fit can be included.
    Errors can be calculated as the sqrt of the histogram counts.
    It returns histogram and fit parameters with their uncertainty.
    It returns also the ax, so as to be able to add more curves to it.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator
    # Calculate the histogram
    h = np.histogram(values, bins=Nbins, range=myrange, density=False)
    y = h[0]
    yerr = np.sqrt(h[0])
    x_edges = h[1]
    x = (x_edges[1:]+x_edges[:-1])*0.5
    dx = x[1] - x[0]
    if IWantDensity:
        y = y / (dx*len(values))
        yerr = yerr / (dx*len(values))
        ylbl = "Density"
    else:
        ylbl = "Counts"
    # Plot
    fig = plt.figure(figsize=myFigSize)
    if icludeErrors and yerr.size>0:
        plt.bar(x, y, width=bw, alpha=0.75) #, label=mylabel)
        plt.errorbar(x, y, yerr, label=mylabel, \
                     linestyle='', linewidth=1, color='b', \
                     marker='o', ms=4, ecolor='b', capsize=0.0)
    else:
        plt.bar(x, y, label=mylabel)
    # Fit
    if IWantGaussianFit:
        def Gauss(x, a, b, c):
            return a*np.exp(-np.power(x - b, 2)/(2*np.power(c, 2)))
        from scipy.optimize import curve_fit
        pars, cov = curve_fit(f=Gauss, xdata=x, ydata=y, \
                              p0=[np.max(y), np.mean(values), np.std(values)], \
                              bounds=(-np.inf, np.inf))
        fit = Gauss(x, *pars)       
        res = y - fit
        stdevs = np.sqrt(np.diag(cov))
        a_value = pars[0]
        mu_value = pars[1]
        sigma_value = pars[2]
        if plotFitPar:
            fitlbl = f'fit: $\\mu$={pars[1]:.2f}+/-{stdevs[1]:.2f}, $\\sigma$={pars[2]:.2f}+/-{stdevs[2]:.2f}'
        else:
            fitlbl = 'fit' 
        plt.plot(x, fit, linestyle='--', linewidth=2, color='black', label=fitlbl)
        for i in range(3):
            print("fit parameter a[%d] of %s: %.2f +/ %.2f" % (i, xlbl, pars[i], stdevs[i]))
        print('\n')
    else:
        pars = []
        stdevs = []
    if (icludeErrors or IWantGaussianFit) and (mylabel != ''):
        plt.legend(fontsize=fs*0.75)
    plt.title(mytitle, fontsize=fs)
    plt.xlabel(xlbl, fontsize=fs)
    plt.ylabel(ylbl, fontsize=fs)
    plt.xticks(fontsize=fs, rotation=0)
    plt.yticks(fontsize=fs, rotation=0)
    plt.gca().xaxis.set_minor_locator(AutoMinorLocator(5))
    plt.gca().yaxis.set_minor_locator(AutoMinorLocator(5))
    plt.gca().tick_params(axis="both", which='major', direction='in', length=8)
    plt.gca().tick_params(axis="both", which='minor', direction='in', length=4)
    #plt.grid(which="major", color="gray", linestyle="--", linewidth=1)
    if saveFig:
        plt.savefig(figname+".jpg")
    ax = plt.gca()
    if showPlot:
        plt.show()
    #Return
    return h, pars, stdevs, ax
    
    
def scatterplot_with_hist(x1, y1, x2=[0], y2=[0], lbl1='', lbl2='', \
                          xlabel='', ylabel='', opacity=1, \
                          nbins=100, myoutpath='', saveFigs=False):
    
    """
    Scatter plot of one or two distributions x1 vs y1 and x2 vs y2. If the distributions
    to plot are actually two, the scatter plots are superimposed.
    It makes also the histograms of x1/x2 and y1/y2 (I can specify the number of bins).
    It returns nothing, but the figure can be saved.
    """
    
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter, MaxNLocator
    from numpy import linspace
    from matplotlib.ticker import AutoMinorLocator
    
    # Set up default x and y limits
    xmm = max([abs(min(x1)), abs(max(x1)), abs(min(x2)), abs(max(x2))])
    ymm = max([abs(min(y1)), abs(max(y1)), abs(min(y2)), abs(max(y2))])
    xlims = [-xmm, xmm]
    ylims = [-ymm, ymm]
    #print("xlims:", xlims)
    #print("ylims:", ylims)

    # Define the locations for the axes
    left, width = 0.12, 0.55
    bottom, height = 0.12, 0.55
    bottom_h = left_h = left + width + 0.02

    # Set up the geometry of the three plots
    rect_main = [left, bottom, width, height] # dimensions of main plot
    rect_histx = [left, bottom_h, width, 0.25] # dimensions of x-histogram
    rect_histy = [left_h, bottom, 0.25, height] # dimensions of y-histogram

    # Set up the size of the figure and other parameters
    fig = plt.figure(1, figsize=(9.5, 9)) 
    fs = 16
    psize = 10
    c1 = '#1f77b4'
    c2 = '#ff7f0e'

    # Make the three plots
    axMain = plt.axes(rect_main) # main plot
    axHistx = plt.axes(rect_histx) # x histogram
    axHisty = plt.axes(rect_histy) # y histogram

    # Make the main plot
    axMain.scatter(x1, y1, s=psize, alpha=opacity, label=lbl1, c=c1)
    if len(x2) > 1: 
        axMain.scatter(x2, y2, s=psize, alpha=opacity, label=lbl2, c=c2)
    axMain.legend(fontsize=fs*0.75)
    
    # Set up the axes/grid style
    axMain.xaxis.set_minor_locator(AutoMinorLocator(5))
    axMain.yaxis.set_minor_locator(AutoMinorLocator(5))
    axMain.tick_params(axis="both", which='major', direction='in', length=8)
    axMain.tick_params(axis="both", which='minor', direction='in', length=4)
    axMain.grid(color='gray', linestyle='--', linewidth=1, alpha=0.5)
    #axMain.axis('equal') #it causes a mismatch between the grids

    # Plot the axes labels
    axMain.set_xlabel(xlabel, fontsize=fs)
    axMain.set_ylabel(ylabel, fontsize=fs)

    # Make the tickmarks pretty
    ticklabels = axMain.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(fs)
        #label.set_family('serif')
    ticklabels = axMain.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(fs)
        #label.set_family('serif')

    # Set up the plot limits
    axMain.set_xlim(xlims)
    axMain.set_ylim(ylims)

    # Set up the histogram bins
    xmin = min(xlims)
    xmax = max(xlims)
    ymin = min(ylims)
    ymax = max(ylims)
    xbins = np.arange(xmin, xmax, (xmax-xmin)/nbins)
    ybins = np.arange(ymin, ymax, (ymax-ymin)/nbins)
    
    # Plot the histograms
    axHistx.hist(x1, bins=xbins, alpha=opacity, label=lbl1, color=c1)
    if len(x2) > 1:
        axHistx.hist(x2, bins=xbins, alpha=opacity, label=lbl2, color=c2)
    #axHistx.legend(fontsize=fs*0.75)
    axHisty.hist(y1, bins=ybins, orientation='horizontal', alpha=opacity, label=lbl1, color=c1)
    if len(x2) > 1:
        axHisty.hist(y2, bins=ybins, orientation='horizontal', alpha=opacity, label=lbl2, color=c2)
    #axHisty.legend(fontsize=fs*0.75)
 
    # Set up the histogram limits
    axHistx.set_xlim(xlims)
    axHisty.set_ylim(ylims)

    # Set up the histogram labels
    axHistx.set_ylabel("Counts [arb. units]", fontsize=fs)
    axHisty.set_xlabel("Counts [arb. units]", fontsize=fs)
    
    # Make the tickmarks pretty
    ticklabels = axHistx.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(fs)
    ticklabels = axHisty.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(fs)

    # Set up the histogram axes/grid style
    axHistx.xaxis.set_minor_locator(AutoMinorLocator(5))
    axHistx.yaxis.set_minor_locator(AutoMinorLocator(5))
    axHistx.tick_params(axis="both", which='major', direction='in', length=8)
    axHistx.tick_params(axis="both", which='minor', direction='in', length=4)
    axHistx.grid(color='gray', linestyle='--', linewidth=1, alpha=0.5)
    axHisty.xaxis.set_minor_locator(AutoMinorLocator(5))
    axHisty.yaxis.set_minor_locator(AutoMinorLocator(5))
    axHisty.tick_params(axis="both", which='major', direction='in', length=8)
    axHisty.tick_params(axis="both", which='minor', direction='in', length=4)
    axHisty.grid(color='gray', linestyle='--', linewidth=1, alpha=0.5)
    axHisty.tick_params(axis='x', labelsize=fs, rotation=-90)
    
    # Remove the inner axes numbers of the histograms
    nullfmt = NullFormatter()
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)    

    # Cool trick that changes the number of tickmarks for the histogram axes
    axHistx.yaxis.set_major_locator(MaxNLocator(4))
    axHisty.xaxis.set_major_locator(MaxNLocator(4))
    
    # Show and save the plot
    if saveFigs:
        plt.savefig(myoutpath + '2D_with_histos' + '.jpg')
    plt.show()    


def plot_hist2d(x, y, \
                NbinX=100, NbinY=100, plot_pixel=True, \
                xlbl='x', ylbl='y', units='mm', mycmap='jet', fs=14, \
                myFigSize=(10, 10), saveFig=False, figname='my2Ddistrib'):
    """
    Function to plot the hist2d of two variables (x,y).
    It reterns the hist2d of the data.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    # Settings
    if plot_pixel:
        xstr = xlbl + ' (%s)' % units
        ystr = ylbl + ' (%s)' % units
        xyrange = [[1, NbinX], [1, NbinY]]
    else:
        xstr = xlbl + ' (pixel)'
        ystr = ylbl + ' (pixel)'
        mmabs = max([abs(min(x)), abs(max(x)), \
                     abs(min(y)), abs(max(y))])
        xyrange = [[-mmabs, mmabs], [-mmabs, mmabs]]
    # Plot
    fig = plt.figure(figsize=myFigSize)
    hh = plt.hist2d(x, y, \
                     range=xyrange, \
                     bins=(NbinX, NbinY), density=False)
    plt.close()
    det_image = hh[0].T
    plt.imshow(det_image, cmap=mycmap, \
               extent=np.array(xyrange).ravel(), \
               interpolation='none')
    cbar = plt.colorbar()    
    cbar.set_label('Counts (arb. units)', fontsize=fs, rotation=90)
    cbar.ax.tick_params(labelsize=fs)
    plt.title("", fontsize=fs)
    plt.xlabel(xstr, fontsize=fs)
    plt.ylabel(ystr, fontsize=fs)
    plt.xticks(fontsize=fs, rotation=0)
    plt.yticks(fontsize=fs, rotation=0)
    plt.axis('equal')
    if saveFig:
        plt.savefig(figname+".jpg")
    plt.show()
    # Return hist2d
    return hh
        

def calc_TH2D_profiles(TH2D, xlimL, xlimH, ylimL, ylimH, XBinEdges, YBinEdges, \
                       p2m=5, plot_mrad=False, IWantFit=False, useP0=False, IWantPlot=False, \
                       lblpltX="data", lblpltY="data", lblfitX="fit", lblfitY="fit", \
                       lblX="", lblY="", lblC="Counts (arb. units)", \
                       plot_title="", lgnd_loc="best", \
                       saveFig=False, outputFile="TH2D_profiles"):

    """
    Function to calculate and plot the profile of a 2D histogram (TH2D).
    It is general, but it was conceived for the deflection distribution
    of a charged particle beam by an oriented crystal. For this reason, 
    the labels set by default refer to this case. If this is not the case,
    provide the proper labels, in particular lblX and lblY.
    It returns the profiles and the parameters of their Gaussian+Linear fit.
    Developed by gpaterno.
    """

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator
    from scipy.optimize import curve_fit
    
    # Sanity check
    if not lblX == "":
        print("Set passed labels! plot_mrad set to False, be sure to pass correct values!\n")
        plot_mrad = False
    
    # Set mult coeff for mrad plots
    if plot_mrad:
        cmr = 1e-3
    else:
        cmr = 1

    # Calculate profiles
    NbinX = len(XBinEdges)-1
    profileX = np.zeros(NbinX)
    ic = int(TH2D.shape[1]/2+1) 
    for j in range(NbinX):
        profileX[j] = np.mean(TH2D[ic-p2m:ic+p2m,j])
    XBinC = XBinEdges[:-1] + (XBinEdges[2]-XBinEdges[1])*0.5

    NbinY = len(YBinEdges)-1
    profileY = np.zeros(NbinY)
    jc = int(TH2D.shape[0]/2+1)
    for i in range(NbinY):
        profileY[i] = np.mean(TH2D[i,jc-p2m:jc+p2m])    
    YBinC = YBinEdges[:-1] + (YBinEdges[2]-YBinEdges[1])*0.5
    
    # Fit profiles
    parsX = []
    parsY = []
    if IWantFit:
        def GaussLin(x, a, b, c, d, e):
            return a*np.exp(-np.power(x - b, 2)/(2*np.power(c, 2))) + (d*x+e)

        try: 
            parsX, covX = curve_fit(f=GaussLin, xdata=XBinC, ydata=profileX, \
                                    #p0=[1, (xlimH+xlimL)*0.5, (xlimH-xlimL)*0.5, 0, 0] if useP0 else None, \
                                    p0=[np.max(profileX), np.mean(profileX), np.std(profileX)] if useP0 else None, \
                                    bounds=(-np.inf, np.inf))
                                    #bounds=(-3*np.std(profileX), 3*np.std(profileX)) if useP0 else (-np.inf, np.inf))    
            fitX = GaussLin(XBinC, *parsX)
            print('Gaussian+Linear fit parsX:\n', parsX)
            resX = profileX - fitX
            stdevsX = np.sqrt(np.diag(covX))
            muX_fit = parsX[1]
            sigmaX_fit = abs(parsX[2])
            
            parsY, covY = curve_fit(f=GaussLin, xdata=YBinC, ydata=profileY, \
                                    #p0=[1, (ylimH+ylimL)*0.5, (ylimH-ylimL)*0.5, 0, 0] if useP0 else None, \
                                    p0=[np.max(profileY), np.mean(profileY), np.std(profileY)] if useP0 else None, \
                                    bounds=(-np.inf, np.inf))
                                    #bounds=(-3*np.std(profileY), 3*np.std(profileY)) if useP0 else (-np.inf, np.inf))   
            fitY = GaussLin(YBinC, *parsY)
            print('Gaussian+Linear fit parsY:\n', parsY, '\n')
            resY = profileY - fitY
            stdevsY = np.sqrt(np.diag(covY))
            muY_fit = parsY[1]
            sigmaY_fit = abs(parsY[2])
        except:
            fitX = 0
            fitY = 0
            sigmaX_fit = 0
            sigmaY_fit = 0
            print('Fit unsuccessful!\n')
    
    # Calculate main statistical values of the selected region
    muX_data, sigmaX_data = weighted_avg_and_std(XBinC, profileX)
    muY_data, sigmaY_data = weighted_avg_and_std(YBinC, profileY)

    if lblX == "":
        if plot_mrad:
            lblpltX = 'data ($\\mu_{data}=%.3f$ mrad, $\\sigma_{data}=%.3f$ mrad)' % \
                      (muX_data*cmr, sigmaX_data*cmr)
            lblpltY = 'data ($\\mu_{data}=%.3f$ mrad, $\\sigma_{data}=%.3f$ mrad)' % \
                      (muY_data*cmr, sigmaY_data*cmr)
            if IWantFit and (sigmaX_fit != 0 and sigmaY_fit != 0):
                lblfitX = 'fit ($\\mu_{fit}=%.3f$ mrad, $\\sigma_{fit}=%.3f$ mrad)' % \
                          (muX_fit*cmr, sigmaX_fit*cmr)
                lblfitY = 'fit ($\\mu_{fit}=%.3f$ mrad, $\\sigma_{fit}=%.3f$ mrad)' % \
                          (muY_fit*cmr, sigmaY_fit*cmr)
            lblX = '$\\Delta\\theta_{x}$ (mrad)'
            lblY = '$\\Delta\\theta_{y}$ (mrad)'
        else:
            lblpltX = 'data ($\\mu_{data}=%.0f$ $\\mu$rad, $\\sigma_{data}=%.0f$ $\\mu$rad)' % \
                      (muX_data*cmr, sigmaX_data*cmr)
            lblpltY = 'data ($\\mu_{data}=%.0f$ $\\mu$rad, $\\sigma_{data}=%.0f$ $\\mu$rad)' % \
                      (muY_data*cmr, sigmaY_data*cmr)
            if IWantFit and (sigmaX_fit != 0 and sigmaY_fit != 0):
                lblfitX = 'fit ($\\mu_{fit}=%.0f$ $\\mu$rad, $\\sigma_{fit}=%.0f$ $\\mu$rad)' % \
                          (muX_fit*cmr, sigmaX_fit*cmr)
                lblfitY = 'fit ($\\mu_{fit}=%.0f$ $\\mu$rad, $\\sigma_{fit}=%.0f$ $\\mu$rad)' % \
                          (muY_fit*cmr, sigmaY_fit*cmr)
            lblX = '$\\Delta\\theta_{x}$ ($\\mu$rad)'
            lblY = '$\\Delta\\theta_{y}$ ($\\mu$rad)'
    else:
        if lblY == "": lblY = "Y"
        lblpltX = 'data ($\\mu_{data}=%.3f$, $\\sigma_{data}=%.3f$)' % (muX_data*cmr, sigmaX_data*cmr)
        lblpltY = 'data ($\\mu_{data}=%.3f$, $\\sigma_{data}=%.3f$)' % (muY_data*cmr, sigmaY_data*cmr)
        if IWantFit and (sigmaX_fit != 0 and sigmaY_fit != 0):
            lblfitX = 'fit ($\\mu_{fit}=%.3f$, $\\sigma_{fit}=%.3f$)' % (muX_fit*cmr, sigmaX_fit*cmr)
            lblfitY = 'fit ($\\mu_{fit}=%.3f$, $\\sigma_{fit}=%.3f$)' % (muY_fit*cmr, sigmaY_fit*cmr)
    
    # Plot profiles
    if IWantPlot:
        fig = plt.figure(figsize=(17, 6))
        fs = 14
        ms = 8
        plt.subplot(1,2,1)
        plt.plot(XBinC*cmr, profileX, linestyle='-', linewidth=2, color='blue', \
                 marker='', markersize=ms, markerfacecolor='blue', \
                 label=lblpltX)
        if IWantFit and (sigmaX_fit != 0 and sigmaY_fit != 0):
            plt.plot(XBinC*cmr, fitX, linestyle='--', linewidth=2, color='black', \
                     label=lblfitX)
        plt.legend(loc=lgnd_loc, fontsize=fs*0.7)
        plt.title(plot_title, fontsize=fs)
        plt.xlabel(lblX, fontsize=fs)
        plt.ylabel(lblC, fontsize=fs, wrap=True)
        plt.xticks(fontsize=fs, rotation=0)
        plt.yticks(fontsize=fs, rotation=0)
        plt.gca().xaxis.set_minor_locator(AutoMinorLocator(5))
        plt.gca().yaxis.set_minor_locator(AutoMinorLocator(5))
        plt.gca().tick_params(axis="both", which='major', direction='in', length=8)
        plt.gca().tick_params(axis="both", which='minor', direction='in', length=4)
        plt.grid('on')
        plt.subplot(1,2,2)
        plt.plot(YBinC*cmr, profileY, linestyle='-', linewidth=2, color='blue', \
                 marker='', markersize=ms, markerfacecolor='blue', \
                 label=lblpltY)
        if IWantFit and (sigmaX_fit != 0 and sigmaY_fit != 0):
            plt.plot(YBinC*cmr, fitY, linestyle='--', linewidth=2, color='black', \
                     label=lblfitY)
        plt.legend(loc=lgnd_loc, fontsize=fs*0.7)
        plt.title(plot_title, fontsize=fs)
        plt.xlabel(lblY, fontsize=fs)
        plt.ylabel(lblC, fontsize=fs, wrap=True)
        plt.xticks(fontsize=fs, rotation=0)
        plt.yticks(fontsize=fs, rotation=0)
        plt.gca().xaxis.set_minor_locator(AutoMinorLocator(5))
        plt.gca().yaxis.set_minor_locator(AutoMinorLocator(5))
        plt.gca().tick_params(axis="both", which='major', direction='in', length=8)
        plt.gca().tick_params(axis="both", which='minor', direction='in', length=4)
        plt.grid('on')
        if saveFig:
            plt.savefig(outputFile + ".jpg")
        plt.show()
        #plt.close()

    # Return the profiles and the fit parameters
    return profileX, parsX, profileY, parsY

    
def find_gaussian_peaks_scipy(x_data, y_data, height=None, distance=None, prominence=None):
    """
    Find peaks using scipy's find_peaks and fit Gaussian curves to them.
    Parameters:
    x_data, y_data: data arrays
    height: minimum height of peaks
    distance: minimum distance between peaks
    prominence: minimum prominence of peaks
    """
    # Import required libraries
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.signal import find_peaks
    from scipy.optimize import curve_fit   
    # Find peaks
    peaks, properties = find_peaks(y_data, height=height, distance=distance, prominence=prominence)
    print(f"Found {len(peaks)} peaks at indices: {peaks}")
    print(f"Peak positions (x): {x_data[peaks]}")
    # Fit Gaussian to each peak
    fitted_params = []
    for i, peak_idx in enumerate(peaks):
        # Define fitting region around peak
        window = min(50, len(x_data)//10)  # Adjust window size as needed
        start = max(0, peak_idx - window)
        end = min(len(x_data), peak_idx + window)
        x_fit = x_data[start:end]
        y_fit = y_data[start:end]
        # Initial guess for parameters
        initial_guess = [y_data[peak_idx], x_data[peak_idx], (x_data[-1]-x_data[0])/20]
        # Do the fit
        try:
            def Gauss(x, amplitude, mean, stddev):
                return amplitude * np.exp(-((x - mean) / stddev)**2 / 2)
            popt, pcov = curve_fit(Gauss, x_fit, y_fit, p0=initial_guess, maxfev=5000)
            fitted_params.append((popt, pcov))
            print(f"Peak {i+1}: amplitude={popt[0]:.3f}, mean={popt[1]:.3f}, stddev={popt[2]:.3f}")
        except Exception as e:
            print(f"Could not fit peak {i+1}: {e}")
            fitted_params.append(None)
    # Return
    return peaks, fitted_params

def create_2d_histogram_and_slice(x_data, y_data, x_bins=50, y_bins=50, XYextent=None, \
                                  x_slice_position=None, slice_width=None, units='', flip_plot=False, \
                                  IWantFit=False, useP0=True, fitTh=100, fitCenter=0, fitExtent=None, \
                                  saveFig=False, figname="TH2D_and_slice"):
    """
    Create a 2D histogram and extract a slice along Y at a given x position with specified width.
    
    Parameters:
    x_data, y_data: Input data arrays.
    x_bins, y_bins: Number of bins for x and y axes.
    XYextent: [range in X, range in Y].
    x_slice_position: X position to take the slice (if None, uses middle of range, which could be != 0).
    slice_width: Width of the slice in x units. If None, uses single bin width.
    units: string with the units of the data (for the plot).
    flip_plot: if True flip the plot w.r.t Y-axis.
    IWantFit: If I want to fit the slice or not.
    useP0: use estimated initial parameters to do the fit.
    fitTh: threshold from which I shift from single Gaussian fit to find peak function.
    fitCenter: suggest the center of the fit.
    fitExtent: range used for the fit.
    
    Returns:
    fig: matplotlib figure object
    y_slice: The slice data along Y (summed over the width)
    y_centers: Y bin centers
    """

    # Import required libraries
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit
    
    # Create 2D histogram
    hist, x_edges, y_edges = np.histogram2d(x_data, y_data, bins=[x_bins, y_bins], range=XYextent)

    # If no slice position specified, use the middle of x range
    if x_slice_position is None:
        x_slice_position = (x_edges[0] + x_edges[-1]) / 2
    
    # Calculate bin width
    x_bin_width = x_edges[1] - x_edges[0]
    
    # If no slice width specified, use single bin width
    if slice_width is None:
        slice_width = x_bin_width
    
    # Calculate how many bins to include in the slice
    n_bins_in_slice = max(1, int(round(slice_width / x_bin_width)))
    
    # Find which x bin contains our slice position
    center_bin_index = np.digitize(x_slice_position, x_edges) - 1
    
    # Calculate start and end bin indices for the slice
    half_width_bins = n_bins_in_slice // 2
    start_bin = center_bin_index - half_width_bins
    end_bin = start_bin + n_bins_in_slice
    
    # Ensure bin indices are within valid range
    start_bin = max(0, start_bin)
    end_bin = min(len(x_edges) - 1, end_bin)
    
    # Extract the slice along Y by summing over the specified x bins
    if start_bin < end_bin:
        y_slice = np.mean(hist[start_bin:end_bin, :], axis=0)
    else:
        # Fallback to single bin if width is too small
        y_slice = hist[center_bin_index, :]
    
    # Calculate actual slice boundaries for plotting
    slice_start = x_edges[start_bin]
    slice_end = x_edges[end_bin]
    actual_width = slice_end - slice_start
    y_centers = (y_edges[:-1] + y_edges[1:]) / 2

    #Define a Gaussian function for fit
    def Gauss(x, a, b, c):
        return a*np.exp(-np.power(x - b, 2)/(2*np.power(c, 2)))

    # Fit profiles
    from scipy.optimize import curve_fit
    x_centers = (x_edges[:-1] + x_edges[1:]) / 2
    ylimH = np.max(x_centers)
    ylimL = np.min(x_centers)
    if not fitExtent:
        p00 = np.max(y_slice)
        p01 = np.mean(y_slice)
        p02 = np.std(y_slice)
        myBounds = (-np.inf, np.inf)
    else:
        useP0 = True
        p00 = np.max(y_slice)
        p01 = fitCenter
        p02 = fitExtent
        myBounds = (fitCenter-fitExtent*0.5, fitCenter+fitExtent*0.5)
        print(myBounds)   
    try:
        if abs(x_slice_position) > fitTh and abs(x_slice_position) < 200:
            sigmaY_fit = 9999
        else:
            parsY, covY = curve_fit(f=Gauss, xdata=y_centers, ydata=y_slice, \
                                    p0=[p00, p01, p02] if useP0 else None, \
                                    bounds=myBounds)    
            fitY = Gauss(y_centers, *parsY)
            resY = y_slice - fitY
            stdevsY = np.sqrt(np.diag(covY))
            muY_fit = parsY[1]
            sigmaY_fit = abs(parsY[2])
    except:
        fitX = 0
        fitY = 0
        sigmaY_fit = 0
        print('Fit unsuccessful!\n')   
    
    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

    if flip_plot:
        matrix_to_plot = np.flip(hist.T, axis=0)
    else:
        matrix_to_plot = hist.T
    
    # Plot 1: 2D histogram
    im = ax1.imshow(matrix_to_plot, origin='lower', aspect='auto', 
                    extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]],
                    cmap='jet')
    
    # Draw the slice region
    ax1.axvline(x=slice_start, color='red', linestyle='--', linewidth=1, alpha=0.7)
    ax1.axvline(x=slice_end, color='red', linestyle='--', linewidth=1, alpha=0.7)
    ax1.axvspan(slice_start, slice_end, alpha=0.3, color='red', \
                label=f'Slice: x={x_slice_position:.2f}±{actual_width/2:.2f}' + units
               )
    ax1.set_aspect('equal')
    ax1.set_xlabel('X ' + '[' + units + ']')
    ax1.set_ylabel('Y ' + '[' + units + ']')
    ax1.set_title('2D Histogram', fontsize=12)
    ax1.legend(fontsize=12)
    plt.colorbar(im, ax=ax1, label='Counts')
    
    # Plot 2: horizontal histogram of the extracted slice
    ax2.barh(y_centers, y_slice, height=y_edges[1]-y_edges[0], alpha=0.7, label='data')
    if IWantFit and sigmaY_fit != 0:
        if abs(x_slice_position) > fitTh and abs(x_slice_position) < 200:
            peaks, fitted_params = \
                find_gaussian_peaks_scipy(y_centers, y_slice, height=0.5, distance=50, prominence=0.3)
            for params in fitted_params:
                if params is not None:
                    popt, _ = params
                    y_fit = gaussian(y_centers, *popt)
                    ax2.plot(y_fit, y_centers, 'k--', alpha=0.7)
        else:
            lblfitY = 'fit ($\\mu_{fit}=%.0f$ %s, $\\sigma_{fit}=%.0f$ %s)' % (muY_fit, units, sigmaY_fit, units)
            ax2.plot(fitY, y_centers, linestyle='--', linewidth=2, color='black', label=lblfitY)       
        ax2.legend(fontsize=10)
    else:
        mean_slice, std_slice = weighted_avg_and_std(y_centers, y_slice)
        lblNoFit = '$\\mu_{data}=%.0f$ %s, $\\sigma_{data}=%.0f$ %s' % (mean_slice, units, std_slice, units)
        xpos_txt = np.max(y_slice)*0.33
        ypos_txt = 2*std_slice        
        plt.text(xpos_txt, ypos_txt, lblNoFit, fontsize=10)
    ax2.set_xlabel('Counts')
    ax2.set_ylabel('Y ' + '[' + units + ']')
    ax2.set_title(f'Slice at x = {x_slice_position:.2f}±{actual_width/2:.2f}' + ' ' + units, fontsize=12)
    ax2.grid(True, alpha=0.3)
    if saveFig:
        plt.savefig(figname + ".jpg")
    plt.tight_layout()
    plt.show()
    
    # Return
    return fig, y_slice, y_centers, (slice_start, slice_end)
#######################################################################################################

############# Set of functions to to save spectra #####################################################
def write_spectrum(Eedges, spectrum, output_file):
    """
    Function to write an energy spectrum obtained through an
    histogram to a text file. Units are arbitrary.
    output_file must be complete: output_path + outname + extension.
    """
    DE = Eedges[1] - Eedges[0]
    Ebin = Eedges[:-1] + DE*0.5
    print("Energy bin width: %.8f" % DE) 
    #Eprob = spectrum * DE
    #print("Spectrum integral:", round(sum(Eprob), 4))
    print("Spectrum integral:", round(sum(spectrum), 4))
    with open(output_file, 'w') as f:
        f.write('#energy spectrum\n')
        for i in range(len(spectrum)):
            #f.write('%.8f %.8f\n' % (Ebin[i], Eprob[i]))
            f.write('%.8f %.8f\n' % (Ebin[i], spectrum[i]))
    print('Spectrum written to %s!\n' % output_file)


def write_spectrum_to_G4file(Eedges_MeV, spectrum, output_file):
    """
    Function to write an energy spectrum obtained through an
    histogram with density=True to a file with Geant4 format.
    In Geant4 the energy bins must be defined in MeV.
    output_file must be complete: output_path + outname + extension.
    """
    DE = Eedges_MeV[1] - Eedges_MeV[0]
    Ebin = Eedges_MeV[:-1] + DE*0.5
    print("Energy bin width: %.8f MeV" % DE)    
    Eprob = spectrum * DE
    print("Spectrum integral (sum(h_i)*DE):", round(sum(Eprob), 4))
    with open(output_file, 'w') as f:
        f.write('#energy spectrum\n')
        f.write('/gps/ene/type User\n')
        f.write('/gps/hist/type energy\n')
        for i in range(len(spectrum)):
            f.write('/gps/hist/point %.8f %.8f\n' % (Ebin[i], Eprob[i]))
    print('Spectrum written to %s!\n' % output_file)
#######################################################################################################    

#######################################################################################################
def get_parula_cmap():
    cm_data_parula = [[0.2422, 0.1504, 0.6603],
    [0.2444, 0.1534, 0.6728],
    [0.2464, 0.1569, 0.6847],
    [0.2484, 0.1607, 0.6961],
    [0.2503, 0.1648, 0.7071],
    [0.2522, 0.1689, 0.7179],
    [0.254, 0.1732, 0.7286],
    [0.2558, 0.1773, 0.7393],
    [0.2576, 0.1814, 0.7501],
    [0.2594, 0.1854, 0.761],
    [0.2611, 0.1893, 0.7719],
    [0.2628, 0.1932, 0.7828],
    [0.2645, 0.1972, 0.7937],
    [0.2661, 0.2011, 0.8043],
    [0.2676, 0.2052, 0.8148],
    [0.2691, 0.2094, 0.8249],
    [0.2704, 0.2138, 0.8346],
    [0.2717, 0.2184, 0.8439],
    [0.2729, 0.2231, 0.8528],
    [0.274, 0.228, 0.8612],
    [0.2749, 0.233, 0.8692],
    [0.2758, 0.2382, 0.8767],
    [0.2766, 0.2435, 0.884],
    [0.2774, 0.2489, 0.8908],
    [0.2781, 0.2543, 0.8973],
    [0.2788, 0.2598, 0.9035],
    [0.2794, 0.2653, 0.9094],
    [0.2798, 0.2708, 0.915],
    [0.2802, 0.2764, 0.9204],
    [0.2806, 0.2819, 0.9255],
    [0.2809, 0.2875, 0.9305],
    [0.2811, 0.293, 0.9352],
    [0.2813, 0.2985, 0.9397],
    [0.2814, 0.304, 0.9441],
    [0.2814, 0.3095, 0.9483],
    [0.2813, 0.315, 0.9524],
    [0.2811, 0.3204, 0.9563],
    [0.2809, 0.3259, 0.96],
    [0.2807, 0.3313, 0.9636],
    [0.2803, 0.3367, 0.967],
    [0.2798, 0.3421, 0.9702],
    [0.2791, 0.3475, 0.9733],
    [0.2784, 0.3529, 0.9763],
    [0.2776, 0.3583, 0.9791],
    [0.2766, 0.3638, 0.9817],
    [0.2754, 0.3693, 0.984],
    [0.2741, 0.3748, 0.9862],
    [0.2726, 0.3804, 0.9881],
    [0.271, 0.386, 0.9898],
    [0.2691, 0.3916, 0.9912],
    [0.267, 0.3973, 0.9924],
    [0.2647, 0.403, 0.9935],
    [0.2621, 0.4088, 0.9946],
    [0.2591, 0.4145, 0.9955],
    [0.2556, 0.4203, 0.9965],
    [0.2517, 0.4261, 0.9974],
    [0.2473, 0.4319, 0.9983],
    [0.2424, 0.4378, 0.9991],
    [0.2369, 0.4437, 0.9996],
    [0.2311, 0.4497, 0.9995],
    [0.225, 0.4559, 0.9985],
    [0.2189, 0.462, 0.9968],
    [0.2128, 0.4682, 0.9948],
    [0.2066, 0.4743, 0.9926],
    [0.2006, 0.4803, 0.9906],
    [0.195, 0.4861, 0.9887],
    [0.1903, 0.4919, 0.9867],
    [0.1869, 0.4975, 0.9844],
    [0.1847, 0.503, 0.9819],
    [0.1831, 0.5084, 0.9793],
    [0.1818, 0.5138, 0.9766],
    [0.1806, 0.5191, 0.9738],
    [0.1795, 0.5244, 0.9709],
    [0.1785, 0.5296, 0.9677],
    [0.1778, 0.5349, 0.9641],
    [0.1773, 0.5401, 0.9602],
    [0.1768, 0.5452, 0.956],
    [0.1764, 0.5504, 0.9516],
    [0.1755, 0.5554, 0.9473],
    [0.174, 0.5605, 0.9432],
    [0.1716, 0.5655, 0.9393],
    [0.1686, 0.5705, 0.9357],
    [0.1649, 0.5755, 0.9323],
    [0.161, 0.5805, 0.9289],
    [0.1573, 0.5854, 0.9254],
    [0.154, 0.5902, 0.9218],
    [0.1513, 0.595, 0.9182],
    [0.1492, 0.5997, 0.9147],
    [0.1475, 0.6043, 0.9113],
    [0.1461, 0.6089, 0.908],
    [0.1446, 0.6135, 0.905],
    [0.1429, 0.618, 0.9022],
    [0.1408, 0.6226, 0.8998],
    [0.1383, 0.6272, 0.8975],
    [0.1354, 0.6317, 0.8953],
    [0.1321, 0.6363, 0.8932],
    [0.1288, 0.6408, 0.891],
    [0.1253, 0.6453, 0.8887],
    [0.1219, 0.6497, 0.8862],
    [0.1185, 0.6541, 0.8834],
    [0.1152, 0.6584, 0.8804],
    [0.1119, 0.6627, 0.877],
    [0.1085, 0.6669, 0.8734],
    [0.1048, 0.671, 0.8695],
    [0.1009, 0.675, 0.8653],
    [0.0964, 0.6789, 0.8609],
    [0.0914, 0.6828, 0.8562],
    [0.0855, 0.6865, 0.8513],
    [0.0789, 0.6902, 0.8462],
    [0.0713, 0.6938, 0.8409],
    [0.0628, 0.6972, 0.8355],
    [0.0535, 0.7006, 0.8299],
    [0.0433, 0.7039, 0.8242],
    [0.0328, 0.7071, 0.8183],
    [0.0234, 0.7103, 0.8124],
    [0.0155, 0.7133, 0.8064],
    [0.0091, 0.7163, 0.8003],
    [0.0046, 0.7192, 0.7941],
    [0.0019, 0.722, 0.7878],
    [0.0009, 0.7248, 0.7815],
    [0.0018, 0.7275, 0.7752],
    [0.0046, 0.7301, 0.7688],
    [0.0094, 0.7327, 0.7623],
    [0.0162, 0.7352, 0.7558],
    [0.0253, 0.7376, 0.7492],
    [0.0369, 0.74, 0.7426],
    [0.0504, 0.7423, 0.7359],
    [0.0638, 0.7446, 0.7292],
    [0.077, 0.7468, 0.7224],
    [0.0899, 0.7489, 0.7156],
    [0.1023, 0.751, 0.7088],
    [0.1141, 0.7531, 0.7019],
    [0.1252, 0.7552, 0.695],
    [0.1354, 0.7572, 0.6881],
    [0.1448, 0.7593, 0.6812],
    [0.1532, 0.7614, 0.6741],
    [0.1609, 0.7635, 0.6671],
    [0.1678, 0.7656, 0.6599],
    [0.1741, 0.7678, 0.6527],
    [0.1799, 0.7699, 0.6454],
    [0.1853, 0.7721, 0.6379],
    [0.1905, 0.7743, 0.6303],
    [0.1954, 0.7765, 0.6225],
    [0.2003, 0.7787, 0.6146],
    [0.2061, 0.7808, 0.6065],
    [0.2118, 0.7828, 0.5983],
    [0.2178, 0.7849, 0.5899],
    [0.2244, 0.7869, 0.5813],
    [0.2318, 0.7887, 0.5725],
    [0.2401, 0.7905, 0.5636],
    [0.2491, 0.7922, 0.5546],
    [0.2589, 0.7937, 0.5454],
    [0.2695, 0.7951, 0.536],
    [0.2809, 0.7964, 0.5266],
    [0.2929, 0.7975, 0.517],
    [0.3052, 0.7985, 0.5074],
    [0.3176, 0.7994, 0.4975],
    [0.3301, 0.8002, 0.4876],
    [0.3424, 0.8009, 0.4774],
    [0.3548, 0.8016, 0.4669],
    [0.3671, 0.8021, 0.4563],
    [0.3795, 0.8026, 0.4454],
    [0.3921, 0.8029, 0.4344],
    [0.405, 0.8031, 0.4233],
    [0.4184, 0.803, 0.4122],
    [0.4322, 0.8028, 0.4013],
    [0.4463, 0.8024, 0.3904],
    [0.4608, 0.8018, 0.3797],
    [0.4753, 0.8011, 0.3691],
    [0.4899, 0.8002, 0.3586],
    [0.5044, 0.7993, 0.348],
    [0.5187, 0.7982, 0.3374],
    [0.5329, 0.797, 0.3267],
    [0.547, 0.7957, 0.3159],
    [0.5609, 0.7943, 0.305],
    [0.5748, 0.7929, 0.2941],
    [0.5886, 0.7913, 0.2833],
    [0.6024, 0.7896, 0.2726],
    [0.6161, 0.7878, 0.2622],
    [0.6297, 0.7859, 0.2521],
    [0.6433, 0.7839, 0.2423],
    [0.6567, 0.7818, 0.2329],
    [0.6701, 0.7796, 0.2239],
    [0.6833, 0.7773, 0.2155],
    [0.6963, 0.775, 0.2075],
    [0.7091, 0.7727, 0.1998],
    [0.7218, 0.7703, 0.1924],
    [0.7344, 0.7679, 0.1852],
    [0.7468, 0.7654, 0.1782],
    [0.759, 0.7629, 0.1717],
    [0.771, 0.7604, 0.1658],
    [0.7829, 0.7579, 0.1608],
    [0.7945, 0.7554, 0.157],
    [0.806, 0.7529, 0.1546],
    [0.8172, 0.7505, 0.1535],
    [0.8281, 0.7481, 0.1536],
    [0.8389, 0.7457, 0.1546],
    [0.8495, 0.7435, 0.1564],
    [0.86, 0.7413, 0.1587],
    [0.8703, 0.7392, 0.1615],
    [0.8804, 0.7372, 0.165],
    [0.8903, 0.7353, 0.1695],
    [0.9, 0.7336, 0.1749],
    [0.9093, 0.7321, 0.1815],
    [0.9184, 0.7308, 0.189],
    [0.9272, 0.7298, 0.1973],
    [0.9357, 0.729, 0.2061],
    [0.944, 0.7285, 0.2151],
    [0.9523, 0.7284, 0.2237],
    [0.9606, 0.7285, 0.2312],
    [0.9689, 0.7292, 0.2373],
    [0.977, 0.7304, 0.2418],
    [0.9842, 0.733, 0.2446],
    [0.99, 0.7365, 0.2429],
    [0.9946, 0.7407, 0.2394],
    [0.9966, 0.7458, 0.2351],
    [0.9971, 0.7513, 0.2309],
    [0.9972, 0.7569, 0.2267],
    [0.9971, 0.7626, 0.2224],
    [0.9969, 0.7683, 0.2181],
    [0.9966, 0.774, 0.2138],
    [0.9962, 0.7798, 0.2095],
    [0.9957, 0.7856, 0.2053],
    [0.9949, 0.7915, 0.2012],
    [0.9938, 0.7974, 0.1974],
    [0.9923, 0.8034, 0.1939],
    [0.9906, 0.8095, 0.1906],
    [0.9885, 0.8156, 0.1875],
    [0.9861, 0.8218, 0.1846],
    [0.9835, 0.828, 0.1817],
    [0.9807, 0.8342, 0.1787],
    [0.9778, 0.8404, 0.1757],
    [0.9748, 0.8467, 0.1726],
    [0.972, 0.8529, 0.1695],
    [0.9694, 0.8591, 0.1665],
    [0.9671, 0.8654, 0.1636],
    [0.9651, 0.8716, 0.1608],
    [0.9634, 0.8778, 0.1582],
    [0.9619, 0.884, 0.1557],
    [0.9608, 0.8902, 0.1532],
    [0.9601, 0.8963, 0.1507],
    [0.9596, 0.9023, 0.148],
    [0.9595, 0.9084, 0.145],
    [0.9597, 0.9143, 0.1418],
    [0.9601, 0.9203, 0.1382],
    [0.9608, 0.9262, 0.1344],
    [0.9618, 0.932, 0.1304],
    [0.9629, 0.9379, 0.1261],
    [0.9642, 0.9437, 0.1216],
    [0.9657, 0.9494, 0.1168],
    [0.9674, 0.9552, 0.1116],
    [0.9692, 0.9609, 0.1061],
    [0.9711, 0.9667, 0.1001],
    [0.973, 0.9724, 0.0938],
    [0.9749, 0.9782, 0.0872],
    [0.9769, 0.9839, 0.0805]]

    from matplotlib.colors import LinearSegmentedColormap
    parula_cmap = LinearSegmentedColormap.from_list('parula', cm_data_parula)
    
    return parula_cmap
#######################################################################################################

#######################################################################################################
def save_dataframe(df, figname, plottitle="", figsize=(1600, 800), header_fontsize=14, cell_fontsize=14):
    """
    Custom wrapper function, based on df2img, to plot a dataframe. 
    """
    import df2img

    fig = df2img.plot_dataframe(
        df.round(2),

        title=dict(
            font_color="black",
            font_family="calibri",
            font_size=24,
            text=plottitle,
            x=0.45,
            xanchor="center",
        ),

        tbl_header=dict(
            align="center",
            fill_color="gray",
            font_color="white",
            font_size=header_fontsize,
            font_family="calibri",
        ),

        tbl_cells=dict(
            height=30,
        ),

        font_size=cell_fontsize,
        row_fill_color=("#ffffff", "#d7d8d6"),
        fig_size=figsize,
    )
    
    df2img.save_dataframe(fig=fig, filename=figname+".jpg")
#######################################################################################################

#######################################################################################################
############# Set of functions to calculate and plot the phase space of particle beams ################
#######################################################################################################
def calc_CAIN_beam_properties(df, m=0.511, beVerbose=True, \
                              IWantSaveStat=False, exportpath='', exportname=''):
    
    """
    Function adapted form GenerateElectronBeamForCain.m (Matlab) script,
    which is based on calculations originally developed by I. Drebot
    for electrons and revised by gpaterno.
    This version accepts a CAIN phase space given in the format:
    columns = ['x[m]', 'y[m]', 't[s]', 'E[eV]', 'px[eV/c]', 'py[eV/c]', 'pz[eV/c]']
    and calculates the properties, returning the norm emittance and the Twiss parameters.
    The emittance calculation is based on (xp, yp), which is valid at small energy spread.
    The particle features can be saved to txt file at a given (passed) path.
    The particle mass in MeV/c2 can be passed. The verbosity can be turned off.
    """
    
    import numpy as np
    import pandas as pd
    
    # constants
    c = 299792458 #speed of light in vacuum [m/s]

    #calculate statistic (trace-based)
    theta = np.arctan((df['px[eV/c]']**2+df['py[eV/c]']**2)**0.5/df['pz[eV/c]']) #[rad]    
    e_x = df['x[m]']                     #[m]
    e_xp = df['px[eV/c]']/df['pz[eV/c]'] #[rad]
    e_y = df['y[m]']                     #[m]
    e_yp = df['py[eV/c]']/df['pz[eV/c]'] #[rad]
    e_s = df['t[s]']*c                   #[m]
    e_E = df['E[eV]']*1e-6               #[MeV]

    e_em_x = np.sqrt(np.mean((e_x-np.mean(e_x))**2)*np.mean((e_xp-np.mean(e_xp))**2)-\
                     np.mean((e_x-np.mean(e_x))*(e_xp-np.mean(e_xp)))**2)
    e_em_y = np.sqrt(np.mean((e_y-np.mean(e_y))**2)*np.mean((e_yp-np.mean(e_yp))**2)-\
                     np.mean((e_y-np.mean(e_y))*(e_yp-np.mean(e_yp)))**2)
    
    mean_x = np.mean(e_x)
    mean_y = np.mean(e_y)
    sigma_x = np.std(e_x)
    sigma_y = np.std(e_y)
    mean_xp = np.mean(e_xp)
    mean_yp = np.mean(e_yp)
    sigma_xp = np.std(e_xp)
    sigma_yp = np.std(e_yp)
    gamma = np.mean(e_E/m)
    delta_gamma = np.std(e_E/m)
    norm_em_x = np.sqrt(gamma**2-1)*e_em_x
    norm_em_y = np.sqrt(gamma**2-1)*e_em_y
    energy_spread = np.std(e_E)/np.mean(e_E)
    mean_e_E = np.mean(e_E)
    std_e_E = np.std(e_E)
    n_m_p = len(e_E)

    #Twiss parameters (gpaterno)
    beta_x = sigma_x**2/e_em_x
    beta_y = sigma_y**2/e_em_y
    gamma_x = sigma_xp**2/e_em_x
    gamma_y = sigma_yp**2/e_em_y
    alpha_x = np.sqrt(beta_x*gamma_x-1)
    alpha_y = np.sqrt(beta_y*gamma_y-1)

    #print results
    if beVerbose:
        print('--------------------------------------------------\n')
        print('Beam non norm Emittances:')
        print('EmitX = %.5e m*rad' % e_em_x)
        print('EmitY = %.5e m*rad' % e_em_y)
        print('Beam normalized Emittances:')
        print('EmitXn = %.5e m*rad' % norm_em_x)
        print('EmitYn = %.5e m*rad' % norm_em_y)
        print('')
        print('sigma_x = %.4f um' % (sigma_x*1e6))
        print('sigma_y = %.4f um' % (sigma_y*1e6))
        print('')
        print('mean_energy = %.4f MeV' % (np.mean(e_E)))
        print('std_energy = %.4f MeV' % (np.std(e_E)))
        print('energy_spread = %.4f' % energy_spread)
        print('gamma = %.4f' % gamma)
        print('delta_gamma = %.4f' % delta_gamma)
        print('--------------------------------------------------\n')
        print('NMP = %d\n' % n_m_p)
        print('theta_mean = %.4e rad' % np.mean(theta))
        print('theta_std = %.4e rad' % np.std(theta))
        print('theta_rms = %.4e rad' % ((np.std(theta)**2+np.mean(theta)**2)**0.5))
        print('ept = %.4e rad' % ((sigma_xp**2+sigma_yp**2)**0.5))
        print('theta_max = %.4e rad' % np.max(theta))
        print('--------------------------------------------------\n')
        print('sigma_z = %.4f um' % np.std(e_s*1e6))
        print('--------------------------------------------------\n')
        print('Twiss parameters:')
        print('alpha_x = %.4f' % alpha_x)
        print('alpha_y = %.4f' % alpha_y)
        print('beta_x = %.6f m' % beta_x)
        print('beta_y = %.6f m' % beta_y)
        print('gamma_x = %.4f m^-1' % gamma_x)
        print('gamma_y = %.4f m^-1' % gamma_y)
        print('--------------------------------------------------\n')

    #export results
    if IWantSaveStat:
        if exportname == '':
            exportname = 'ebeam_stat.txt'
        with open(exportpath+exportname, 'w') as f:
            f.write('--------------------------------------------------\n')
            f.write('Beam non norm Emittances:\n')
            f.write('EmitX = %.5e m*rad\n' % e_em_x)
            f.write('EmitY = %.5e m*rad\n' % e_em_y)
            f.write('Beam normalized Emittances:\n')
            f.write('EmitXn = %.5e m*rad\n' % norm_em_x)
            f.write('EmitYn = %.5e m*rad\n' % norm_em_y)
            f.write('\n')
            f.write('sigma_x = %.4f um\n' % (sigma_x*1e6))
            f.write('sigma_y = %.4f um\n' % (sigma_y*1e6))
            f.write('\n')
            f.write('mean_energy = %.4f MeV\n' % (np.mean(e_E)))
            f.write('std_energy = %.4f MeV\n' % (np.std(e_E)))
            f.write('energy_spread = %.4f\n' % energy_spread)
            f.write('gamma = %.4f\n' % gamma)
            f.write('delta_gamma = %.4f\n' % delta_gamma)
            f.write('--------------------------------------------------\n')
            f.write('NMP = %d\n' % n_m_p)
            f.write('theta_mean = %.4e rad\n' % np.mean(theta))
            f.write('theta_std = %.4e rad\n' % np.std(theta))
            f.write('theta_rms = %.4e rad\n' % ((np.std(theta)**2+np.mean(theta)**2)**0.5))
            f.write('ept = %.4e rad\n' % ((sigma_xp**2+sigma_yp**2)**0.5))
            f.write('theta_max = %.4e rad\n' % np.max(theta))
            f.write('--------------------------------------------------\n')
            f.write('sigma_z = %.4f um\n' % np.std(e_s*1e6))
            f.write('--------------------------------------------------\n')
            f.write('Twiss parameters:\n')
            f.write('alpha_x = %.4f\n' % alpha_x)
            f.write('alpha_y = %.4f\n' % alpha_y)
            f.write('beta_x = %.6f m\n' % beta_x)
            f.write('beta_y = %.6f m\n' % beta_y)
            f.write('gamma_x = %.4f m^-1\n' % gamma_x)
            f.write('gamma_y = %.4f m^-1\n' % gamma_y)
            f.write('--------------------------------------------------\n')
        print('Beam statistics written in %s file!\n' % (exportpath+exportname))
    
    #return normalized emittance and Twiss parameters
    return [norm_em_x, norm_em_y], [alpha_x, alpha_y, beta_x, beta_y, gamma_x, gamma_y]


def calc_RFTrack_beam_properties(df, m=0.511, beVerbose=True, \
                                 IWantSaveStat=False, exportpath='', exportname=''):
    
    """
    Function adapted form GenerateElectronBeamForCain.m (Matlab) script,
    which is based on calculations originally developed by I. Drebot
    for electrons and revised by gpaterno.
    This version accepts an RF-Track phase space given in the format:
    columns = ["x[mm]", "xp[mrad]", "y[mm]", "yp[mrad]", "t[mm/c]", "p[MeV/c]", "ID"]
    and calculates the properties, returning the norm emittance and the Twiss parameters.
    The emittance calculation is based on (px, py), which is the most general (RF-Track).
    The particle features can be saved to txt file at a given (passed) path.
    The particle mass in MeV/c2 can be passed. Verbosity can be turned off.
    """
    
    import numpy as np
        
    #calculate statistics    
    x = df['x[mm]'].values*1e-3           #m
    xp = df['xp[mrad]'].values*1e-3       #rad
    y = df['y[mm]'].values*1e-3           #m
    yp = df['yp[mrad]'].values*1e-3       #rad
    bl = df['t[mm/c]'].values             #mm
    p = df['p[MeV/c]'].values             #MeV/c
    E = np.sqrt(p**2+m**2)                #MeV
    pz = p / np.sqrt(xp**2+yp**2+1)       #MeV/c
    px = xp*pz                            #MeV/c
    py = yp*pz                            #MeV/c
    theta = np.arctan((xp**2+yp**2)**0.5) #rad
    gamma = np.mean(E/m)
    delta_gamma = np.std(E/m)
    
    em_x_tr = np.sqrt(np.mean((x-np.mean(x))**2)*np.mean((xp-np.mean(xp))**2)-\
                      np.mean((x-np.mean(x))*(xp-np.mean(xp)))**2)
    em_y_tr = np.sqrt(np.mean((y-np.mean(y))**2)*np.mean((yp-np.mean(yp))**2)-\
                      np.mean((y-np.mean(y))*(yp-np.mean(yp)))**2)
    norm_em_x_tr = np.sqrt(gamma**2-1)*em_x_tr
    norm_em_y_tr = np.sqrt(gamma**2-1)*em_y_tr
    # Previous (trace-based) definition holds only when the beam energy spread is small.
    # A more geneal definition is given by Floettmann PRAB, 034202 (2003).
    # It is used by RF-Track and ASTRA, so I adopted it.
    em_x = np.sqrt(np.mean((x-np.mean(x))**2)*np.mean((px-np.mean(px))**2)-\
                   np.mean((x-np.mean(x))*(px-np.mean(px)))**2)/(m*np.mean(pz))
    em_y = np.sqrt(np.mean((y-np.mean(y))**2)*np.mean((py-np.mean(py))**2)-\
                   np.mean((y-np.mean(y))*(py-np.mean(py)))**2)/(m*np.mean(pz))
    norm_em_x = np.mean(pz)*em_x
    norm_em_y = np.mean(pz)*em_y
    
    mean_x = np.mean(x)
    mean_y = np.mean(y)
    sigma_x = np.std(x)
    sigma_y = np.std(y)
    mean_xp = np.mean(xp)
    mean_yp = np.mean(yp)
    sigma_xp = np.std(xp)
    sigma_yp = np.std(yp)
    mean_px = np.mean(px)
    mean_py = np.mean(py)
    sigma_px = np.std(px)
    sigma_py = np.std(py)
    mean_E = np.mean(E)
    sigma_E = np.std(E)
    energy_spread = sigma_E/mean_E
    n_m_p = len(E)

    #Twiss parameters (gpaterno)
    beta_x = sigma_x**2/em_x_tr
    beta_y = sigma_y**2/em_y_tr
    gamma_x = sigma_xp**2/em_x_tr
    gamma_y = sigma_yp**2/em_y_tr
    xpx_rms = np.mean((x-np.mean(x))*(xp-np.mean(xp)))
    ypy_rms = np.mean((y-np.mean(y))*(yp-np.mean(yp)))
    alpha_x = -xpx_rms/em_x_tr
    alpha_y = -ypy_rms/em_y_tr
    #alpha_x = np.sqrt(beta_x*gamma_x-1) #it is equivalent, if we use em_x_tr
    #alpha_y = np.sqrt(beta_y*gamma_y-1) #it is equivalent, if we use em_y_tr

    #print results
    if beVerbose:
        print('--------------------------------------------------\n')
        print('Beam trace-space Emittances:')
        print('em_x_tr = %.2e m*rad' % (em_x_tr))
        print('em_y_tr = %.2e m*rad' % (em_y_tr))
        print('Beam trace-space Emittances:')
        print('em_x_n_tr = %.2e m*rad' % (norm_em_x_tr))
        print('em_y_n_tr = %.2e m*rad' % (norm_em_y_tr))
        print('Beam Emittances:')
        print('em_x = %.2e m*rad' % (em_x))
        print('em_y = %.2e m*rad' % (em_y))
        print('Beam normalized Emittances:')
        print('em_x_n = %.2e m*rad' % (norm_em_x))
        print('em_y_n = %.2e m*rad' % (norm_em_y))
        print('')
        print('sigma_x = %.2f mm' % (sigma_x*1e3))
        print('sigma_y = %.2f mm' % (sigma_y*1e3))
        print('sigma_xp = %.2f mrad' % (sigma_xp*1e3))
        print('sigma_yp = %.2f mrad' % (sigma_yp*1e3))
        print('sigma_px = %.2f MeV/c' % sigma_px)
        print('sigma_py = %.2f MeV/c' % sigma_py)
        print('')
        print('theta_mean = %.4e rad' % np.mean(theta))
        print('theta_std = %.4e rad' % np.std(theta))
        print('theta_rms = %.4e rad' % ((np.std(theta)**2+np.mean(theta)**2)**0.5))
        print('theta_max = %.4e rad' % np.max(theta))
        print('--------------------------------------------------\n')
        print('bunch_length = %.2f mm' % np.std(bl))
        print('')
        print('mean_energy = %.4f MeV' % mean_E)
        print('std_energy = %.4f MeV' % sigma_E)
        print('energy_spread = %.4f' % energy_spread)
        print('gamma = %.4f' % gamma)
        print('delta_gamma = %.4f' % delta_gamma)
        print('--------------------------------------------------\n')
        print('NMP = %d\n' % n_m_p)
        print('--------------------------------------------------\n')
        print('Twiss parameters (calculated using em_x,y_tr):')
        print('alpha_x = %.4f' % alpha_x)
        print('alpha_y = %.4f' % alpha_y)
        print('beta_x = %.6f m' % beta_x)
        print('beta_y = %.6f m' % beta_y)
        print('gamma_x = %.4f m^-1' % gamma_x)
        print('gamma_y = %.4f m^-1' % gamma_y)
        print('--------------------------------------------------\n')
        
    #export results
    if IWantSaveStat:
        if exportname == '':
            exportname = 'beam_stat.txt'
        with open(exportpath+exportname, 'w') as f:
            f.write('--------------------------------------------------\n')
            f.write('Beam trace-space Emittances:')
            f.write('em_x_tr = %.2e m*rad' % (em_x_tr))
            f.write('em_y_tr = %.2e m*rad' % (em_y_tr))
            f.write('Beam trace-space Emittances:')
            f.write('em_x_n_tr = %.2e m*rad' % (norm_em_x_tr))
            f.write('em_y_n_tr = %.2e m*rad' % (norm_em_y_tr))
            f.write('Beam Emittances:')
            f.write('em_x = %.2e m*rad' % (em_x))
            f.write('em_y = %.2e m*rad' % (em_y))
            f.write('Beam normalized Emittances:')
            f.write('em_x_n = %.2e m*rad' % (norm_em_x))
            f.write('em_y_n = %.2e m*rad' % (norm_em_y))
            f.write('')
            f.write('sigma_x = %.2f mm' % (sigma_x*1e3))
            f.write('sigma_y = %.2f mm' % (sigma_y*1e3))
            f.write('sigma_xp = %.2f mrad' % (sigma_xp*1e3))
            f.write('sigma_yp = %.2f mrad' % (sigma_yp*1e3))
            f.write('sigma_px = %.2f MeV/c' % sigma_px)
            f.write('sigma_py = %.2f MeV/c' % sigma_py)
            f.write('')
            f.write('theta_mean = %.4e rad' % np.mean(theta))
            f.write('theta_std = %.4e rad' % np.std(theta))
            f.write('theta_rms = %.4e rad' % ((np.std(theta)**2+np.mean(theta)**2)**0.5))
            f.write('theta_max = %.4e rad' % np.max(theta))
            f.write('--------------------------------------------------\n')
            f.write('bunch_length = %.2f mm' % np.std(bl))
            f.write('')
            f.write('mean_energy = %.4f MeV' % mean_E)
            f.write('std_energy = %.4f MeV' % sigma_E)
            f.write('energy_spread = %.4f' % energy_spread)
            f.write('gamma = %.4f' % gamma)
            f.write('delta_gamma = %.4f' % delta_gamma)
            f.write('--------------------------------------------------\n')
            f.write('NMP = %d\n' % n_m_p)
            f.write('--------------------------------------------------\n')
            f.write('Twiss parameters (calculated using em_x,y_tr):')
            f.write('alpha_x = %.4f' % alpha_x)
            f.write('alpha_y = %.4f' % alpha_y)
            f.write('beta_x = %.6f m' % beta_x)
            f.write('beta_y = %.6f m' % beta_y)
            f.write('gamma_x = %.4f m^-1' % gamma_x)
            f.write('gamma_y = %.4f m^-1' % gamma_y)
            f.write('--------------------------------------------------\n')            
        print('Beam statistics written in %s file!\n' % (exportpath+exportname))
    
    #return normalized emittance and Twiss parameters
    return [norm_em_x, norm_em_y], [alpha_x, alpha_y, beta_x, beta_y, gamma_x, gamma_y]


def plot_RFTrack_transverse_phase_space(df, m=0.511, radius_sel=1e15, \
                                        num_bins_x=50, num_bins_xp=50, \
                                        num_bins_y=50, num_bins_yp=50, \
                                        xrange=None, xprange=None, \
                                        yrange=None, yprange=None, \
                                        use_log_scale=False, mymap='jet', \
                                        myoutpath='', saveFigs=False):
    
    """
    It accepts a dataframe with an RF-Track phase space given in the format:
    columns = ["x[mm]", "xp[mrad]", "y[mm]", "yp[mrad]", "t[mm/c]", "p[MeV/c]", "ID"]
    and plots the 2d histogram of "x[mm]" and "xp[mrad], and the same for y.
    Also, the spatial distribution of beam intensity and energy is plotted. 
    A circular selection can be applied to the dataframe.
    """    
    
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors

    df_sel = df[df['x[mm]']**2 + df['y[mm]']**2 < radius_sel**2] #apply a selection
 
    plt.figure(figsize=(18, 6))
    fs = 16
    plt.subplot(1,2,1)
    hist_matrix1, x_edges, xp_edges = np.histogram2d(df_sel['x[mm]'], df_sel['xp[mrad]'], \
                                                     bins=(num_bins_x, num_bins_xp), \
                                                     range=[xrange, xprange]
                                                    )  
    img1 = plt.imshow(hist_matrix1.T, cmap=mymap, interpolation="None", aspect='auto', origin='lower', \
                      extent=[x_edges[0], x_edges[-1], xp_edges[0], xp_edges[-1]]
                     )
    if use_log_scale:
        img1.set_norm(mcolors.LogNorm())
    cbar = plt.colorbar(img1, ax=plt.gca(), label='Frequency')
    cbar.set_label('counts (arb. units)', fontsize=fs, rotation=90)
    cbar.ax.tick_params(labelsize=fs)
    plt.xlabel('x (mm)', fontsize=fs)
    plt.ylabel('xp (mrad)', fontsize=fs)
    plt.xticks(fontsize=fs, rotation=0)
    plt.yticks(fontsize=fs, rotation=0)
    plt.subplot(1,2,2)
    hist_matrix2, y_edges, yp_edges = np.histogram2d(df_sel['y[mm]'], df_sel['yp[mrad]'], \
                                                     bins=(num_bins_y, num_bins_yp), \
                                                     range=[yrange, yprange]
                                                    )  
    img2 = plt.imshow(hist_matrix2.T, cmap=mymap, interpolation="None", aspect='auto', origin='lower', 
                      extent=[y_edges[0], y_edges[-1], yp_edges[0], yp_edges[-1]]
                     )
    if use_log_scale:
        img2.set_norm(mcolors.LogNorm())
    cbar = plt.colorbar(img2, ax=plt.gca(), label='Frequency')
    cbar.set_label('counts (arb. units)', fontsize=fs, rotation=90)
    cbar.ax.tick_params(labelsize=fs)
    plt.xlabel('y (mm)', fontsize=fs)
    plt.ylabel('yp (mrad)', fontsize=fs)
    plt.xticks(fontsize=fs, rotation=0)
    plt.yticks(fontsize=fs, rotation=0)
    if saveFigs:
        plt.savefig(myoutpath + 'xVSxpANDyVSyp' + '.jpg')
    plt.show()  
    
    xp_sel = df_sel['xp[mrad]'].values
    yp_sel = df_sel['yp[mrad]'].values
    E_sel = np.sqrt(df_sel['p[MeV/c]']**2 + m**2).values
    #df_sel['E[MeV]'] = np.sqrt(df_sel['p[MeV/c]']**2 + m**2)    
    Zmatrix, Zmatrix_norm, XbinEdges, YbinEdges = TProfile2D(xp_sel, yp_sel, \
                                                             num_bins_xp, *xprange, \
                                                             num_bins_yp, *yprange, \
                                                             E_sel)
    
    plt.figure(figsize=(18, 6))
    fs = 16
    plt.subplot(1,2,1)
    hist_matrix3, _, _ = np.histogram2d(df_sel['x[mm]'], df_sel['y[mm]'], \
                                        bins=(num_bins_x, num_bins_y), \
                                        range=[xrange, yrange]
                                       )  
    img3 = plt.imshow(hist_matrix3.T, cmap=mymap, interpolation="None", aspect='auto', origin='lower', \
                      extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]]
                     )
    if use_log_scale:
        img3.set_norm(mcolors.LogNorm())
    cbar = plt.colorbar(img3, ax=plt.gca(), label='Frequency')
    cbar.set_label('counts (arb. units)', fontsize=fs, rotation=90)
    cbar.ax.tick_params(labelsize=fs)
    plt.xlabel('x (mm)', fontsize=fs)
    plt.ylabel('y (mm)', fontsize=fs)
    plt.xticks(fontsize=fs, rotation=0)
    plt.yticks(fontsize=fs, rotation=0)
    plt.axis('equal')
    plt.subplot(1,2,2)
    Nlev = 50
    #plt.contourf(XbinEdges, YbinEdges, Zmatrix_norm, Nlev, cmap=mymap)
    plt.imshow(Zmatrix_norm, extent=[xp_edges[0], xp_edges[-1], yp_edges[0], yp_edges[-1]], cmap=mymap)
    cbar = plt.colorbar()
    cbar.set_label('E (MeV)', fontsize=fs, rotation=90)
    cbar.ax.tick_params(labelsize=fs)
    plt.xlabel("$\\theta_{x}$ (mrad)", fontsize=fs)
    plt.ylabel("$\\theta_{y}$ (mrad)", fontsize=fs)
    plt.xticks(fontsize=fs, rotation=0)
    plt.yticks(fontsize=fs, rotation=0)
    plt.axis('equal')    
    if saveFigs:
        plt.savefig(myoutpath + 'IandE' + '.jpg')
    plt.show()  
    
    
def plot_RFTrack_longitudinal_phase_space(df, \
                                          num_bins_x=50, num_bins_y=50, myrange=None, \
                                          use_log_scale=False, mymap='jet', \
                                          myoutpath='', saveFigs=False):

    """
    It accepts a dataframe with an RF-Track phase space given in the format:
    columns = ["x[mm]", "xp[mrad]", "y[mm]", "yp[mrad]", "t[mm/c]", "p[MeV/c]", "ID"]
    and plots the 2d histogram of "t[mm/c]" and "p[MeV/c].
    """    
    
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors

    fig, ax = plt.subplots(figsize=(14, 10))
    fs = 18
    hist_matrix, x_edges, y_edges = np.histogram2d(df['t[mm/c]'], df['p[MeV/c]'], \
                                                   bins=(num_bins_x, num_bins_y), \
                                                   range=myrange, \
                                                  )
    img = ax.imshow(hist_matrix.T, cmap=mymap, interpolation="None", aspect='auto', origin='lower', \
                    extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]]
                   )
    if use_log_scale:
        img.set_norm(mcolors.LogNorm())
    cbar = plt.colorbar(img, ax=ax, label='Frequency')
    cbar.set_label('counts (arb. units)', fontsize=fs, rotation=90)
    cbar.ax.tick_params(labelsize=fs)
    ax.set_xlabel('t (mm/c)', fontsize=fs)
    ax.set_ylabel('p (MeV/c)', fontsize=fs)
    plt.xticks(fontsize=fs, rotation=45)
    plt.yticks(fontsize=fs, rotation=0)
    if saveFigs: 
        plt.savefig(myoutpath + 'pVSz' + '.jpg')
    plt.show()

    
def plot_RFTrack_xpx_ypy(df, m=0.511, radius_sel=1e15, \
                         angleLim=200, use_log_scale=False, \
                         myoutpath='', saveFigs=False):

    """
    It accepts a dataframe with an RF-Track phase space given in the format:
    columns = ["x[mm]", "xp[mrad]", "y[mm]", "yp[mrad]", "t[mm/c]", "p[MeV/c]", "ID"]
    and do the scatter plot of "x[mm]" and "xp[mrad] (and the same for y)
    weighted by energy. A circular selection can be applied to the dataframe.
    """ 
    
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    df_sel = df[df['x[mm]']**2 + df['y[mm]']**2 < radius_sel**2] #apply a selection
    
    if not use_log_scale:
        colours = df_sel['p[MeV/c]']
        clbl = r"p (MeV/c)"
    else:
        colours = np.log(df_sel['p[MeV/c]'])
        clbl = r"log(p[MeV/c])"
    
    plt.figure(figsize=(18, 5))
    fs = 16
    psize = 10
    opacity = 0.4
    plt.subplot(1,2,1)
    plt.scatter(df_sel['x[mm]'], df_sel['xp[mrad]'], s=psize, alpha=opacity, c=colours)
    cbar = plt.colorbar(ax=plt.gca())
    cbar.ax.tick_params(labelsize=fs)
    cbar.set_label(clbl, fontsize=fs, rotation=90)
    plt.title("", fontsize=fs)
    plt.xlabel("x (mm)", fontsize=fs)
    plt.ylabel("xp (mrad)", fontsize=fs)
    plt.xticks(fontsize=fs, rotation=0)
    plt.yticks(fontsize=fs, rotation=0)
    plt.ylim(-angleLim, angleLim)
    plt.subplot(1,2,2)
    plt.scatter(df_sel['y[mm]'], df_sel['yp[mrad]'], s=psize, alpha=opacity, c=colours)
    cbar = plt.colorbar(ax=plt.gca())
    cbar.ax.tick_params(labelsize=fs)
    cbar.set_label(clbl, fontsize=fs, rotation=90)
    plt.title("", fontsize=fs)
    plt.xlabel("y (mm)", fontsize=fs)
    plt.ylabel("yp (mrad)", fontsize=fs)
    plt.xticks(fontsize=fs, rotation=0)
    plt.yticks(fontsize=fs, rotation=0)
    plt.ylim(-angleLim, angleLim)
    if saveFigs:
        plt.savefig(myoutpath + 'xpx_ypy_w' + '.jpg')
    plt.show()


def plot_RFTrack_x_vs_z(df, lbl, weight):
    
    """
    It accepts a dataframe with an RF-Track phase space given in the format:
    columns = ["x[mm]", "xp[mrad]", "y[mm]", "yp[mrad]", "t[mm/c]", "p[MeV/c]", "ID"]
    and do the scatter plot of "x[mm]" vs "t[mm/c] weighted by another variable of the
    dataframe (e.g. xp[mrad] or p[MeV/c]).
    """ 
    
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    
    plt.figure(figsize=(18, 8))
    fs = 16
    psize = 10
    opacity = 0.4
    plt.scatter(df['t[mm/c]'], df['x[mm]'], label=lbl, s=psize, alpha=opacity, c=np.log(df[weight]))
    cbar = plt.colorbar(ax=plt.gca())
    cbar.ax.tick_params(labelsize=fs)
    cbar.set_label(r"log(%s)" % weight, fontsize=fs, rotation=90)
    plt.legend(fontsize=fs)
    plt.title("", fontsize=fs)
    plt.xlabel("t (mm/c)", fontsize=fs)
    plt.ylabel("x (mm)", fontsize=fs)
    plt.xticks(fontsize=fs, rotation=0)
    plt.yticks(fontsize=fs, rotation=0)
    #plt.axis('equal')
    plt.ylim(-80, 80)
    plt.show()   


def plot_RFTrack_weighted_scatterplots(df1, df2, radius_sel=1e15, lbl1='', lbl2='', \
                                       weightThem=False, weightVar='log(p)', \
                                       myoutpath='', saveFigs=False):

    """
    It accepts two dataframes with an RF-Track phase space given in the format:
    columns = ["x[mm]", "xp[mrad]", "y[mm]", "yp[mrad]", "t[mm/c]", "p[MeV/c]", "ID"]
    and do the scatter plot of "x1[mm]" vs "y1[mm] and "x2[mm]" vs "y2[mm]", eventually
    weighted by the natural logarithm of the particle momentum "p[MeV/c]".
    A circular selection can be applied to the dataframes.
    """ 
    
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    #apply a selection (particle within an aperture)
    if not df1.empty:
        df1_sel = df1[df1['x[mm]']**2 + df1['y[mm]']**2 < radius_sel**2] 
    if not df2.empty:
        df2_sel = df2[df2['x[mm]']**2 + df2['y[mm]']**2 < radius_sel**2]
        
    if weightThem:
        if weightVar == 'log(p)':
            c1 = np.log(df1_sel['p[MeV/c]'])
            clbl = r"log(p[MeV/c])"
        else:
            c1 = (df1_sel['t[mm/c]'])
            clbl = r"t[mm/c]"
        if not df2.empty:
            if weightVar == 'log(p)':
                c2 = np.log(df2_sel['p[MeV/c]'])
            else:
                c2 = (df2_sel['t[mm/c]'])
    else:
        c1 = '#1f77b4'
        c2 = '#ff7f0e'
    
    #plot
    plt.figure(figsize=(8, 8))
    fs = 16
    psize = 10
    opacity = 0.4
    if not df1.empty:
        plt.scatter(df1_sel['x[mm]'], df1_sel['y[mm]'], s=psize, alpha=opacity, label=lbl1, c=c1)
    if not df2.empty:
        plt.scatter(df2_sel['x[mm]'], df2_sel['y[mm]'], s=psize, alpha=opacity, label=lbl2, c=c2)
    if weightThem:
        cbar = plt.colorbar(ax=plt.gca())
        cbar.set_label(clbl, fontsize=fs, rotation=90)
    plt.legend(fontsize=fs)
    plt.title("", fontsize=fs)
    plt.xlabel("x (mm)", fontsize=fs)
    plt.ylabel("y (mm)", fontsize=fs)
    plt.xticks(fontsize=fs, rotation=0)
    plt.yticks(fontsize=fs, rotation=0)
    plt.axis('equal')
    if saveFigs:
        plt.savefig(myoutpath + 'xy' + '.jpg')
    plt.show()    


def plot_RFTrack_transverse_distributions(df, m=0.511, radius_sel=1e15, \
                                          num_bins_x=50, num_bins_xp=50, \
                                          num_bins_y=50, num_bins_yp=50, \
                                          xrange=None, xprange=None, \
                                          yrange=None, yprange=None, \
                                          use_log_scale=False, mymap='jet', \
                                          ttlSD='', ttlE='', \
                                          myoutpath='', saveFigs=False):

    df_sel = df[df['x[mm]']**2 + df['y[mm]']**2 < radius_sel**2] #apply a selection

    xp_sel = df_sel['xp[mrad]'].values
    yp_sel = df_sel['yp[mrad]'].values
    E_sel = np.sqrt(df_sel['p[MeV/c]']**2 + m**2).values
    Zmatrix, Zmatrix_norm, XbinEdges, YbinEdges = TProfile2D(xp_sel, yp_sel, \
                                                             num_bins_xp, *xprange, \
                                                             num_bins_yp, *yprange, \
                                                             E_sel)

    plt.figure(figsize=(18, 6))
    fs = 16
    plt.subplot(1,2,1)
    hist_matrix3, x_edges, y_edges = np.histogram2d(df_sel['x[mm]'], df_sel['y[mm]'], \
                                        bins=(num_bins_x, num_bins_y), \
                                        range=[xrange, yrange])  
    img3 = plt.imshow(hist_matrix3.T, cmap=mymap, interpolation="None", aspect='auto', origin='lower', \
                      extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]])
    if use_log_scale:
        img3.set_norm(mcolors.LogNorm())
    cbar = plt.colorbar(img3, ax=plt.gca(), label='Frequency')
    cbar.set_label('counts (arb. units)', fontsize=fs, rotation=90)
    cbar.ax.tick_params(labelsize=fs)
    plt.title(ttlSD, fontsize=fs)
    plt.xlabel('x (mm)', fontsize=fs)
    plt.ylabel('y (mm)', fontsize=fs)
    plt.xticks(fontsize=fs, rotation=0)
    plt.yticks(fontsize=fs, rotation=0)
    plt.axis('equal')
    plt.subplot(1,2,2)
    Nlev = 50
    plt.contourf(XbinEdges, YbinEdges, Zmatrix_norm, Nlev, cmap=mymap)
    cbar = plt.colorbar()
    cbar.set_label('E (MeV)', fontsize=fs, rotation=90)
    cbar.ax.tick_params(labelsize=fs)
    plt.title(ttlE, fontsize=fs)
    plt.xlabel("$\\theta_{x}$ (mrad)", fontsize=fs)
    plt.ylabel("$\\theta_{y}$ (mrad)", fontsize=fs)
    plt.xticks(fontsize=fs, rotation=0)
    plt.yticks(fontsize=fs, rotation=0)
    plt.axis('equal')    
    if saveFigs:
        plt.savefig(myoutpath + 'IandE' + '.jpg')
    plt.show()

    return hist_matrix3.T, Zmatrix_norm


def plot_EorPspectrum(p, lbl='', \
                      p2=[], lbl2='', \
                      isE=False, plotLog=False, opacity=1, \
                      IWantDensity=False, NormMax=False, \
                      nbin_E=100, range_E=None, solidPlot=False, \
                      myoutpath='', saveFigs=False):
    
    """
    It takes one or two vectors of p [MeV/c] and plot the histogram.
    It can plot momentum p or energy E by setting properly isE variable.
    """

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator
    
    if isE:     
        xlbl = 'E [MeV]'
    else:
        xlbl = 'p [MeV/c]'
        
    if IWantDensity:
        if isE:     
            ylbl = '1/N$\\times$dN/dE (1/MeV)'
        else:
            ylbl = '1/N$\\times$dN/dp (1/MeV/c)'        
        NormMax = False
    else:
        ylbl = 'Counts (arb. units)'
               
    fig = plt.figure(figsize=(9, 6))
    fs = 16
    lw = 2
    
    if solidPlot:
        spectrum, edges, _ = plt.hist(p, density=IWantDensity, bins=nbin_E, \
                                      range=range_E, alpha=opacity, label=lbl)
        if len(p2)>0:     
            spectrum2, _, _ = plt.hist(p2, density=IWantDensity, bins=nbin_E, \
                                       range=range_E, alpha=opacity, label=lbl2)       
    else:
        # code to plot only lines and not solid histograms 
        spectrum, edges = np.histogram(p, density=IWantDensity, bins=nbin_E, range=range_E)
        bin_E = edges[:-1] + (edges[1]-edges[0])*0.5
        if len(p2)>0:
            spectrum2, _ = np.histogram(p2, density=IWantDensity, bins=nbin_E, range=range_E)

        if NormMax:
            spectrum = spectrum / np.max(spectrum)
            if len(p2)>0:
                spectrum2 = spectrum2 / np.max(spectrum2)
        
        plt.plot(bin_E, spectrum, linewidth=lw, alpha=opacity, label=lbl)
        if len(p2)>0:
            plt.plot(bin_E, spectrum2, linewidth=lw, alpha=opacity, label=lbl2)
        
    plt.legend(fontsize=fs*0.75)
    plt.xlabel(xlbl, fontsize=fs)
    plt.ylabel(ylbl, fontsize=fs)
    plt.xticks(fontsize=fs, rotation=0)
    plt.yticks(fontsize=fs, rotation=0)
    plt.gca().xaxis.set_minor_locator(AutoMinorLocator(5))
    plt.gca().yaxis.set_minor_locator(AutoMinorLocator(5))
    plt.gca().tick_params(axis="both", which='major', direction='in', length=8)
    plt.gca().tick_params(axis="both", which='minor', direction='in', length=4)
    #plt.xlim([0, 1000])
    if plotLog:
        plt.yscale('log')
    plt.grid(which="major", color="gray", linestyle="--", linewidth=1)
    plt.title('')
    if saveFigs:
        plt.savefig(myoutpath + 'spectrum' + '.jpg')
    plt.show()
    
    
def plot_EvsTheta(theta, energy, \
                  num_bins_theta=50, num_bins_energy=50, \
                  theta_range=None, energy_range=None,\
                  use_log_scale=False, mymap='jet', \
                  xlbl='$\\theta$ (mrad)', ttl='', \
                  myoutpath='', saveFigs=False):
    
    """
    It takes theta/thetax [mrad] and energy in [MeV] and plot a 2Dhist,
    namely a mustage plot in the case of ICS.
    """

    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    
    fig, ax = plt.subplots(figsize=(8, 5))
    fs = 16
    hist_matrix, x_edges, y_edges = np.histogram2d(theta, energy, \
                                                   bins=(num_bins_theta, num_bins_energy), \
                                                   range=(theta_range, energy_range), \
                                                  )
    img = ax.imshow(hist_matrix.T, cmap=mymap, interpolation="None", aspect='auto', origin='lower', 
                    extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]]
                   )
    if use_log_scale:
        img.set_norm(mcolors.LogNorm())
    cbar = plt.colorbar(img, ax=ax, label='Frequency')
    cbar.set_label('counts (arb. units)', fontsize=fs, rotation=90)
    cbar.ax.tick_params(labelsize=fs)
    ax.set_title(ttl, fontsize=fs)
    ax.set_xlabel(xlbl, fontsize=fs)
    ax.set_ylabel('E (MeV)', fontsize=fs)
    plt.xticks(fontsize=fs, rotation=0)
    plt.yticks(fontsize=fs, rotation=0)
    #ax.invert_yaxis()
    if saveFigs: 
        plt.savefig(myoutpath + 'EvsTheta' + '.jpg')
    plt.show()
#######################################################################################################

