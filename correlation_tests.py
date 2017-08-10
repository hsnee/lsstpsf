import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns;sns.set_style('darkgrid')

x = np.linspace(-3,3,1001)
y = np.linspace(-3,3,1001)
X, Y = np.meshgrid(x, y)
ra = X.flatten()
dec = Y.flatten()
g1 = ra
g2 = dec

# Now instead of plotting things let's input them into treecorr and get a gamma-gamma correlation function.
import treecorr

cat = treecorr.Catalog(g1=g1, g2=g2, ra=ra, dec=dec,ra_units='radians',dec_units='radians')
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
plt.legend([r'$\xi_+$',r'$\xi_-$'],bbox_to_anchor=(0.98, 1), loc='upper right', borderaxespad=0.)
plt.savefig('diverging_test_higher.pdf')
plt.close()

for i in range(len(ra)):
    print ra[i],dec[i],g1[i],g2[i]
