from angles import r2d
import sys, os, matplotlib, lsst.daf.persistence,  galsim
import numpy as np
from lsst.afw.geom import Point2D

if len(sys.argv)<2:
    sys.stderr.write("Syntax: python test-psf.py  repo_path\n")
    exit()

repo_rel_path = sys.argv[1]
repo_abs_path = os.path.abspath(repo_rel_path)

if not os.path.exists(repo_abs_path):
    sys.stderr.write("Nothing found at {}\n".format(repo_abs_path))
    exit()

matplotlib.use("pdf")
import matplotlib.pyplot as plt
import seaborn as sns;sns.set_style('darkgrid')

# Create a data butler which provides access to a data repository.
butler = lsst.daf.persistence.Butler(repo_abs_path)
print "Butler summoned (i.e. we have loaded the data repository)", repo_abs_path
ccd_exposures = butler.queryMetadata("calexp", format=['visit', 'ccd'])
ccd_exposures = [(visit,ccd) for (visit,ccd) in ccd_exposures if butler.datasetExists("calexp", visit=visit, ccd=ccd)]

print "Found {} (exposure,ccd) pairs".format(len(ccd_exposures))
for ccd in range(104):
    if ccd ==1:
        old_X, old_Y, old_g1, old_g2 = X, Y, g1, g2
    visit, ccd = ccd_exposures[ccd]
    calexp = butler.get("calexp", visit=visit, ccd=ccd, immediate=True)
    detector = calexp.getDetector()
    print 'visit', visit, 'ccd', ccd, detector.getName()
    width, height, psf, wcs = calexp.getWidth(), calexp.getHeight(), calexp.getPsf(), calexp.getWcs()
    nx = 20
    ny = int(nx * (height/width))
    X,Y,g1,g2 = (np.zeros((ny,nx)) for _ in range(4))
    for i,x in enumerate(np.linspace(0, width, nx)):
        for j,y in enumerate(np.linspace(0, height, ny)):
            point = Point2D(x,y)
            image = psf.computeKernelImage(point)
            galsim_image = galsim.Image(image.getArray())
            shape_data = galsim.hsm.FindAdaptiveMom(galsim_image,strict=False)
            if shape_data.error_message=='':
                pass
            else:
                print shape_data.error_message
                continue
            FpPos = wcs.pixelToSky(point)
            X[j,i], Y[j,i] = FpPos.getLongitude(), FpPos.getLatitude()
            g1[j,i], g2[j,i] = shape_data.observed_shape.g1, shape_data.observed_shape.g2
    if ccd == 0:
        pass
    else:
        old_X, old_Y, old_g1, old_g2 = np.append(X,old_X), np.append(Y,old_Y), np.append(g1,old_g1), np.append(g2,old_g2)

# Now instead of plotting things let's input them into treecorr and get a gamma-gamma correlation function.
import treecorr
del X,Y,g1,g2
ra = old_X.copy()
dec = old_Y.copy()
g1 = old_g1.copy()
g2 = old_g2.copy()
cat = treecorr.Catalog(g1=g1, g2=g2, ra=ra, dec=dec,ra_units='radians',dec_units='radians')
gg = treecorr.GGCorrelation(min_sep=1, max_sep=200, nbins=40, sep_units='arcmin')
gg.process(cat)
xip = gg.xip
xim = gg.xim

# plot the correlation functions:
plt.figure()
r = np.exp(gg.meanlogr)
plt.xlabel(r'$\theta$ (arcmin)')
plt.ylabel(r'$\xi$')
plt.plot(r,xip,label=r'$\xi_+$')
plt.plot(r,xim,label=r'$\xi_-$')
plt.xscale('log')
plt.yscale('log')
plt.legend(bbox_to_anchor=(0.9, 1), loc='upper right', borderaxespad=0.)
plt.savefig('GGCorrelation.pdf')
plt.close()
