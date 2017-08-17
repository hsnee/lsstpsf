import numpy as np
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
import seaborn as sns;sns.set_style('darkgrid')

ra  = np.random.randn(int(1E5))*12000 - 6000
dec = np.random.randn(int(1E5))*12000 - 6000

theta = np.arctan(dec/ra)
r = np.sqrt(ra**2+dec**2)
g1 = r*np.cos(2*theta)
g2 = r*np.sin(2*theta)

# Now instead of plotting things let's input them into treecorr and get a gamma-gamma correlation function.
import treecorr

cat = treecorr.Catalog(g1=g1, g2=g2, ra=ra1, dec=dec1,ra_units='radians',dec_units='radians')
gg = treecorr.GGCorrelation(min_sep=1, max_sep=200, nbins=40, sep_units='arcmin')
gg.process(cat)
xip = gg.xip
xim = gg.xim
sigma = gg.varxi**0.5

# plot the correlation functions:
r = np.exp(gg.meanlogr)
