from angles import r2arcs, r2d
import sys, os, matplotlib
from lsst.afw.geom import Point2D

if len(sys.argv)<3:
    sys.stderr.write("Syntax: python test-psf.py  repo_path visitnum\n")
    exit()

repo_rel_path = sys.argv[1]
repo_abs_path = os.path.abspath(repo_rel_path)

if not os.path.exists(repo_abs_path):
    sys.stderr.write("Nothing found at {}\n".format(repo_abs_path))
    exit()
visitnum = sys.argv[2]
visitnum = int(visitnum)

import numpy as np
matplotlib.use("pdf")
import matplotlib.pyplot as plt
import lsst.daf.persistence,  galsim
import seaborn as sns;sns.set_style('darkgrid')

# Create a data butler which provides access to a data repository.
butler = lsst.daf.persistence.Butler(repo_abs_path)
print "Butler summoned (i.e. we have loaded the data repository)", repo_abs_path
ccd_exposures = butler.queryMetadata("calexp", format=['visit', 'ccd'])
ccd_exposures = [(visit,ccd) for (visit,ccd) in ccd_exposures if butler.datasetExists("calexp", visit=visit, ccd=ccd)]

print "Found {} (exposure,ccd) pairs".format(len(ccd_exposures))
for ccd in range(visitnum*112,visitnum*112+104):
    if ccd ==1+112*visitnum:
        old_X, old_Y, old_U, old_C, old_angles = X, Y, U, C, angles
    visit, ccd = ccd_exposures[ccd]
    print 'visit', visit, 'ccd', ccd
    calexp = butler.get("calexp", visit=visit, ccd=ccd, immediate=True)
    detector = calexp.getDetector()
    print detector.getName()
    width, height, psf, wcs = calexp.getWidth(), calexp.getHeight(), calexp.getPsf(), calexp.getWcs()

    # pixel_scale = wcs.pixelScale() # if using pixel coordinates.
    pixel_scale = 1. # if using FOCAL_PLANE coordinates
    nx = 5
    ny = int(nx * (height/width))
    angles = np.zeros((ny,nx))
    X,Y,U,C = (np.zeros((ny,nx)) for _ in range(4))
    for i,x in enumerate(np.linspace(0, width, nx)):
        for j,y in enumerate(np.linspace(0, height, ny)):
            point = Point2D(x,y)
            image = psf.computeKernelImage(point)
            array = image.getArray()
            galsim_image = galsim.Image(array)
            shape_data = galsim.hsm.FindAdaptiveMom(galsim_image,strict=False)
            if shape_data.error_message=='':
                pass
            else:
                print shape_data.error_message
                continue
            FpPos = wcs.pixelToSky(point)
            X[j,i], Y[j,i] = FpPos.getLongitude(), FpPos.getLatitude()
            e1, e2 = shape_data.observed_shape.e1, shape_data.observed_shape.e2
            U[j,i] = np.sqrt(e1**2+e2**2)
            C[j,i] = shape_data.moments_sigma
            if C[j,i]<0:
                print 'negative size detected :/'
            angles[j,i] = r2d(0.5*np.arctan2(e2,e1))
    if ccd == 0:
        pass
    else:
        old_X, old_Y, old_U, old_C, old_angles = np.append(X,old_X), np.append(Y,old_Y), np.append(U,old_U), np.append(C,old_C), np.append(angles, old_angles)

old_V = np.zeros(np.shape(old_U))
old_X1, old_Y1, old_U1, old_V1, old_C1, old_angles1 = old_X.copy(), old_Y.copy(), old_U.copy(), old_V.copy(), old_C.copy(), old_angles.copy()
del X,Y,U,V,C,angles,old_X,old_Y,old_U,old_V,old_C,old_angles
# now do the same for the second visit
visitnum=1
for ccd in range(visitnum*112,visitnum*112+104):
    if ccd ==1+112*visitnum:
        old_X, old_Y, old_U, old_C, old_angles = X, Y, U, C, angles
    visit, ccd = ccd_exposures[ccd]
    print 'visit', visit, 'ccd', ccd
    calexp = butler.get("calexp", visit=visit, ccd=ccd, immediate=True)
    detector = calexp.getDetector()
    print detector.getName()
    width, height, psf, wcs = calexp.getWidth(), calexp.getHeight(), calexp.getPsf(), calexp.getWcs()

    # pixel_scale = wcs.pixelScale() # if using pixel coordinates.
    pixel_scale = 1. # if using FOCAL_PLANE coordinates
    nx = 5
    ny = int(nx * (height/width))
    angles = np.zeros((ny,nx))
    X,Y,U,C = (np.zeros((ny,nx)) for _ in range(4))
    for i,x in enumerate(np.linspace(0, width, nx)):
        for j,y in enumerate(np.linspace(0, height, ny)):
            point = Point2D(x,y)
            image = psf.computeKernelImage(point)
            array = image.getArray()
            galsim_image = galsim.Image(array)
            shape_data = galsim.hsm.FindAdaptiveMom(galsim_image,strict=False)
            if shape_data.error_message=='':
                pass
            else:
                print shape_data.error_message
                continue
            FpPos = wcs.pixelToSky(point)
            X[j,i], Y[j,i] = FpPos.getLongitude(), FpPos.getLatitude()
            e1, e2 = shape_data.observed_shape.e1, shape_data.observed_shape.e2
            U[j,i] = np.sqrt(e1**2+e2**2)
            C[j,i] = shape_data.moments_sigma
            if C[j,i]<0:
                print 'negative size detected :/'
            angles[j,i] = r2d(0.5*np.arctan2(e2,e1))
    if ccd == 0:
        pass
    else:
        old_X, old_Y, old_U, old_C, old_angles = np.append(X,old_X), np.append(Y,old_Y), np.append(U,old_U), np.append(C,old_C), np.append(angles, old_angles)

old_V = np.zeros(np.shape(old_U))
old_X2, old_Y2, old_U2, old_V2, old_C2, old_angles2 = old_X.copy(), old_Y.copy(), old_U.copy(), old_V.copy(), old_C.copy(), old_angles.copy()


std2fwhm = 2.*np.sqrt(2.*np.log(2.))
pixel_scale = r2arcs(pixel_scale)
plt.figure()
Q = plt.quiver(pixel_scale*(old_X1-old_X2),pixel_scale*(old_Y1-old_Y2),old_U1-old_U2,old_V1-old_V2,std2fwhm*(old_C1-old_C2,angles=old_angles1-old_angles2,headlength=0,headaxislength=0,scale=20,cmap='viridis')
qk = plt.quiverkey(Q, 0.9, 0.9, 0.05, r'$|e|=0.05$', labelpos='E', coordinates='figure')
plt.title("whisker plot")
plt.colorbar(label=r'$\Delta$ FWHM')
#plt.axis('equal')
plt.xlabel('arcsec')
plt.ylabel('arcsec')
plt.savefig("delta-whisker.pdf")
plt.close()
