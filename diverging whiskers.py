# really diverging

X = np.random.randn(int(1E5))
Y = np.random.randn(int(1E5))
r = np.sqrt(X**2+Y**2)
theta = np.arctan(Y/X)
e1 = r*np.cos(2*theta)
e2 = r*np.sin(2*theta)
U = np.sqrt(e1**2+e2**2)
angles = r2d(0.5*np.arctan2(e2,e1))

# plotting
std2fwhm = 2.*np.sqrt(2.*np.log(2.))
V = np.zeros(np.shape(U))
pixel_scale = r2arcs(1)

plt.figure()
Q = plt.quiver(pixel_scale*X,pixel_scale*Y,U,V,std2fwhm*C,angles=angles,headlength=0,headaxislength=0,scale=20,cmap='viridis')
qk = plt.quiverkey(Q, 0.9, 0.9, 0.05, r'$|e|=0.05$', labelpos='E', coordinates='figure')
plt.title("whisker plot")
plt.colorbar(label='FWHM')
plt.legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0.)
plt.xlabel('arcsec')
plt.ylabel('arcsec')
plt.savefig("really_diverging_whisker.pdf")
plt.close()

# not really diverging

e1 = X
e2 = Y
U = np.sqrt(e1**2+e2**2)
angles = r2d(0.5*np.arctan2(e2,e1))

#plotting
std2fwhm = 2.*np.sqrt(2.*np.log(2.))
V = np.zeros(np.shape(U))
pixel_scale = r2arcs(1)

plt.figure()
Q = plt.quiver(pixel_scale*X,pixel_scale*Y,U,V,std2fwhm*C,angles=angles,headlength=0,headaxislength=0,scale=20,cmap='viridis')
qk = plt.quiverkey(Q, 0.9, 0.9, 0.05, r'$|e|=0.05$', labelpos='E', coordinates='figure')
plt.title("whisker plot")
plt.colorbar(label='FWHM')
plt.legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0.)
plt.xlabel('arcsec')
plt.ylabel('arcsec')
plt.savefig("not_really_diverging_whisker.pdf")
plt.close()
