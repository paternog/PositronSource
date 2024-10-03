# graphical settings for fancier plots

import mplhep as hep

style = hep.style.ROOT

# figure size
style['figure.figsize'] = (10.795, 7.000)

# other sizes
style['lines.markersize'] = 6
style['lines.linewidth'] = 3

# fonts
style['font.family'] = "serif"
style['axes.titlesize'] = "x-small"
style['axes.labelsize'] = "small"
style['xtick.labelsize'] = "x-small"
style['ytick.labelsize'] = "x-small"
style['legend.fontsize'] = "x-small"
style['mathtext.fontset'] = "dejavuserif"

# text locations
style['xaxis.labellocation'] = 'center'
style['yaxis.labellocation'] = 'center'

# grid
#style["axes.grid"] = True

# output plots
out_format = "pdf"
dpi = 800