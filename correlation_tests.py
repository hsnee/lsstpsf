import numpy as np
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
import seaborn as sns;sns.set_style('darkgrid')

data0 = np.loadtxt('0 corr_arrays.txt')
data1 = np.loadtxt('1 corr_arrays.txt')
data2 = np.loadtxt('2 corr_arrays.txt')
ra0,dec0,g10,g20 = data0[0],data0[1],data0[2],data0[3]
ra1,dec1,g11,g21 = data1[0],data1[1],data1[2],data1[3]
ra2,dec2,g12,g22 = data2[0],data2[1],data2[2],data2[3]
g1_avg = 0.5*(g10 + g11)
g2_avg = 0.5*(g20 + g21)
delta_g1 = g11-g10
delta_g2 = g21-g20
# Now instead of plotting things let's input them into treecorr and get a gamma-gamma correlation function.
import treecorr

cat = treecorr.Catalog(g1=delta_g1/g1_avg, g2=delta_g2/g2_avg, ra=ra1, dec=dec1,ra_units='radians',dec_units='radians')
gg = treecorr.GGCorrelation(min_sep=1, max_sep=200, nbins=40, sep_units='arcmin')
gg.process(cat)
xip = gg.xip
xim = gg.xim
sigma = gg.varxi**0.5

# plot the correlation functions:
import seaborn.timeseries

def _plot_std_bars(std=None, central_data=None, ci=None, data=None,*args, **kwargs):
    std = sigma
    ci = np.asarray((central_data - std, central_data + std))
    kwargs.update({"central_data": central_data, "ci": ci, "data": data})
    seaborn.timeseries._plot_ci_band(*args, **kwargs)
seaborn.timeseries._plot_std_bars = _plot_std_bars
r = np.exp(gg.meanlogr)
plt.figure()
sns.tsplot(xip,r,err_style='std_bars')
sns.tsplot(xim,r,err_style='std_bars',color='r')
plt.xlabel(r'$\theta$ (arcmin)')
plt.ylabel(r'$\xi$')
plt.xscale('log')
plt.yscale('log')
plt.legend([r'$\xi_+$',r'$\xi_-$'],bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0.)
plt.savefig('fractional_delta_shear.pdf')
plt.close()
