import numpy as np
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
import lsst.daf.persistence,  galsim
import seaborn as sns;sns.set_style('darkgrid')
from angles import r2arcs, r2d

# really diverging

X = np.random.randn(int(1E4))
Y = np.random.randn(int(1E4))
r = np.sqrt(X**2+Y**2)
theta = np.arctan(Y/X)
angles = np.zeros(len(X))
e1 = r*np.cos(2*theta)
e2 = r*np.sin(2*theta)
U = np.sqrt(e1**2+e2**2)
for i in range(len(e1)):
    angles[i] = r2d(0.5*np.arctan2(e2[i],e1[i]))

# plotting
std2fwhm = 2.*np.sqrt(2.*np.log(2.))
V = np.zeros(np.shape(U))
pixel_scale = r2arcs(1)

plt.figure()
Q = plt.quiver(pixel_scale*X,pixel_scale*Y,U,V,angles=angles,headlength=0,headaxislength=0,scale=20,cmap='viridis')
qk = plt.quiverkey(Q, 0.9, 0.9, 0.05, r'$|e|=0.05$', labelpos='E', coordinates='figure')
plt.xlabel('arcsec')
plt.ylabel('arcsec')
plt.savefig("really_diverging_whisker.pdf")
plt.close()

# not really diverging
del U,V,e1,e2,angles
e1 = X.copy()
e2 = Y.copy()
U = np.sqrt(e1**2+e2**2)
angles = np.zeros(len(X))
for i in range(len(e1)):
    angles[i] = r2d(0.5*np.arctan2(e2[i],e1[i]))
#plotting
std2fwhm = 2.*np.sqrt(2.*np.log(2.))
V = np.zeros(np.shape(U))
pixel_scale = r2arcs(1)

plt.figure()
Q = plt.quiver(pixel_scale*X,pixel_scale*Y,U,V,angles=angles,headlength=0,headaxislength=0,scale=20,cmap='viridis')
qk = plt.quiverkey(Q, 0.9, 0.9, 0.05, r'$|e|=0.05$', labelpos='E', coordinates='figure')
plt.xlabel('arcsec')
plt.ylabel('arcsec')
plt.savefig("not_really_diverging_whisker.pdf")
plt.close()
