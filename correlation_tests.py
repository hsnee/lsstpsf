def arcmin2arcsec(x):
    return x*60
import numpy as np
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
import seaborn as sns;sns.set_style('darkgrid')
from plotting import HSC_style_plots as HSC_plt
from plotting import whisker_plot, sns_time_series

data0 = np.loadtxt('0 corr_arrays.txt')
data1 = np.loadtxt('1 corr_arrays.txt')


ra0,dec0,g10,g20,sigma0 = np.array(map(lambda x: data0[x], (0,1,2,3,4)))
ra1,dec1,g11,g21,sigma1 = np.array(map(lambda x: data1[x], (0,1,2,3,4)))


if 0:
    # get the g-g correlation function from treecorr
    import treecorr
    cat = treecorr.Catalog(g1=delta_g1/g1_avg, g2=delta_g2/g2_avg, ra=ra1, dec=dec1,ra_units='radians',dec_units='radians')
    gg = treecorr.GGCorrelation(min_sep=1, max_sep=200, nbins=40, sep_units='arcmin')
    gg.process(cat)
    xip = gg.xip
    xim = gg.xim
    sigma = gg.varxi**0.5

if 1:
    # call a function from the plotting module
    HSC_plt(ra0,dec0,g10,g20,g11,g21,sigma0,sigma1)
