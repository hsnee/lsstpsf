from angles import r2arcs
import sys, os, matplotlib
import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.geom as geom

if len(sys.argv)<2:
    sys.stderr.write("Syntax: python test-psf.py  repo_path\n")
    sys.exit(1)

repo_rel_path = sys.argv[1]
repo_abs_path = os.path.abspath(repo_rel_path)

if not os.path.exists(repo_abs_path):
    sys.stderr.write("Nothing found at {}\n".format(repo_abs_path))
    sys.exit(1)


import numpy as np
matplotlib.use("pdf")
import matplotlib.pyplot as plt
import lsst.daf.persistence,  galsim
import seaborn as sns;sns.set_style('darkgrid')

# Create a data butler which provides access to a data repository.
butler = lsst.daf.persistence.Butler(repo_abs_path)
print "Butler summoned (i.e. we have loaded the data repository)"
print repo_abs_path
ccd_exposures = butler.queryMetadata("calexp", format=['visit', 'ccd'])
ccd_exposures = [(visit,ccd) for (visit,ccd) in ccd_exposures if butler.datasetExists("calexp", visit=visit, ccd=ccd)]

print "Found {} (exposure,ccd) pairs".format(len(ccd_exposures))
for ccd in range(104):
    if ccd ==1:
        old_X = X
        old_Y = Y
        old_U = U
        old_C = C
        old_angles = angles
    else:
        pass
    if ccd==0:
        pass
    else:
        old_pixel_scale = pixel_scale
    visit, ccd = ccd_exposures[ccd]
    print 'visit', visit, 'ccd', ccd
    calexp = butler.get("calexp", visit=visit, ccd=ccd, immediate=True)
    detector = calexp.getDetector()
    print detector.getName()
    width = calexp.getWidth()
    height = calexp.getHeight()
    psf = calexp.getPsf()
    wcs = calexp.getWcs()
    #pixel_scale = wcs.pixelScale() # if using pixel coordinates.
    pixel_scale = 1. # if using FOCAL_PLANE coordinates
    if ccd==0:
        pass
    else:
        np.testing.assert_approx_equal(old_pixel_scale, pixel_scale, significant=5, err_msg='pixel scale is different!!')
    nx = 5
    ny = int(nx * (height/width))
    angles = np.zeros((ny,nx))
    X,Y,U,C = (np.zeros((ny,nx)) for _ in range(4))
    for i,x in enumerate(np.linspace(0, width, nx)):
        for j,y in enumerate(np.linspace(0, height, ny)):
            image = psf.computeKernelImage(lsst.afw.geom.Point2D(x, y))
            array = image.getArray()
            galsim_image = galsim.Image(array)
            shape_data = galsim.hsm.FindAdaptiveMom(galsim_image,strict=False)
            if shape_data.error_message=='':
                pass
            else:
                print shape_data.error_message
                continue
            FpPos = wcs.pixelToSky(lsst.afw.geom.Point2D(x,y))
            X[j,i] = FpPos.getLongitude()
            Y[j,i] = FpPos.getLatitude()
            e1 = shape_data.observed_shape.e1
            e2 = shape_data.observed_shape.e2
            U[j,i] = np.sqrt(e1**2+e2**2)
            C[j,i] = shape_data.moments_sigma
            if C[j,i]<0:
                print 'negative size detected :/'
            else:
                pass
            angles[j,i] = 0.5*np.arctan2(e2,e1)*180/np.pi
    if ccd == 0:
        pass
    else:
        old_X = np.append(X,old_X)
        old_Y = np.append(Y,old_Y)
        old_U = np.append(U,old_U)
        old_C = np.append(C,old_C)
        old_angles = np.append(angles, old_angles)

# Plotting
std2fwhm = 2.*np.sqrt(2.*np.log(2.))
old_V = np.zeros(np.shape(old_U))
pixel_scale = r2arcs(pixel_scale)

plt.figure()
Q = plt.quiver(pixel_scale*(old_X-np.median(old_X)),pixel_scale*(old_Y-np.median(old_Y)),old_U,old_V,std2fwhm*old_C,angles=old_angles,headlength=0,headaxislength=0,scale=20,cmap='viridis')
qk = plt.quiverkey(Q, 0.9, 0.9, 0.05, r'$|e|=0.05$', labelpos='E',
                   coordinates='figure')
plt.title("whisker plot")
plt.colorbar(label='FWHM')
#plt.axis('equal')
plt.xlabel('arcsec')
plt.ylabel('arcsec')
plt.savefig("/global/homes/h/husni/lsstpsf/whisker.pdf")
plt.close()
