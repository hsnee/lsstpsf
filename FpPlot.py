#import seaborn as sns;sns.set_style('darkgrid')
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
matplotlib.use("agg")
import matplotlib.pyplot as plt
import lsst.daf.persistence,  galsim


# Create a data butler which provides access to a data repository.
butler = lsst.daf.persistence.Butler(repo_abs_path)
print "Butler summoned (i.e. we have loaded the data repository)"
print repo_abs_path
ccd_exposures = butler.queryMetadata("calexp", format=['visit', 'ccd'])
ccd_exposures = [(visit,ccd) for (visit,ccd) in ccd_exposures if butler.datasetExists("calexp", visit=visit, ccd=ccd)]

print "Found {} (exposure,ccd) pairs".format(len(ccd_exposures))
print 0
#nx = 100
#ny = 200
#angles = np.zeros((112,ny,nx))
#X = np.zeros((112,ny,nx))
#Y = np.zeros((112,ny,nx))
#U = np.zeros((112,ny,nx))
#V = np.zeros((112,ny,nx))
#C = np.zeros((112,ny,nx))
pixOx = [[0 for _ in range(2)] for _ in range(112)]
j= 0
#print 0
#positions = np.zeros(112)
#print 0
previousvisit=0
pixOx = np.array(pixOx)
#print 0
savelast=0
for ccd in [6,107]:
    visit, ccd = ccd_exposures[ccd]
    print 'visit', visit, 'ccd', ccd
    calexp = butler.get("calexp", visit=visit, ccd=ccd, immediate=True)
    detector = calexp.getDetector()
    orientation = detector.getOrientation()
    print orientation.getFpPosition()
    print detector.getName()
    pixOx = detector.getCenter(cameraGeom.FOCAL_PLANE)
    pixOx = pixOx.getPoint()
    #print pixOx                                                                                                                                                            
    #print orientation.getReferencePoint()                                                                                                                                  
    #print orientation.getNQuarter(),orientation.getYaw(),orientation.getPitch(),orientation.getRoll()                                                                      
    width = calexp.getWidth()
    height = calexp.getHeight()
    psf = calexp.getPsf()
    wcs = calexp.getWcs()
    pixel_scale = wcs.pixelScale()
    help(psf)
for ccd in range(104):
    #print ccd
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
    orientation = detector.getOrientation()
    print orientation.getFpPosition()
    print detector.getName()
    pixOx = detector.getCenter(cameraGeom.FOCAL_PLANE)
    pixOx = pixOx.getPoint()
    #print pixOx
    #print orientation.getReferencePoint()
    #print orientation.getNQuarter(),orientation.getYaw(),orientation.getPitch(),orientation.getRoll()
    width = calexp.getWidth()
    height = calexp.getHeight()
    psf = calexp.getPsf()
    wcs = calexp.getWcs()
    #pixel_scale = wcs.pixelScale()
    pixel_scale = 1.
    if ccd==0:
        pass
    else: 
        if old_pixel_scale == pixel_scale:
            pass
        else:
            print 'pixel scale is different!!!',old_pixel_scale, pixel_scale
    nx = 5
    ny = int(nx * (height/width))
    angles = np.zeros((ny,nx))
    X = np.zeros((ny,nx))
    Y = np.zeros((ny,nx))
    #print pixOx
    U = np.zeros((ny,nx))
    C = np.zeros((ny,nx))
    for i,x in enumerate(np.linspace(0, width, nx)):
        for j,y in enumerate(np.linspace(0, height, ny)):
            image = psf.computeKernelImage(lsst.afw.geom.Point2D(x, y))
            array = image.getArray()
            galsim_image = galsim.Image(array)
            if ccd==107:
                plt.figure()
                plt.imshow(array)
                plt.savefig(str(ccd)+'.png')
                plt.close()
            else:
                pass
            shape_data = galsim.hsm.FindAdaptiveMom(galsim_image,strict=False)
            if shape_data.error_message=='':
                pass
            else:
                print shape_data.error_message
                continue
            extent = geom.Extent2D(x,y)
            pixO = wcs.pixelToSky(lsst.afw.geom.Point2D(x,y))
            #print pixO,'0000'
            #help(pixOx)
            X[j,i] = pixO.getLongitude()
            Y[j,i] = pixO.getLatitude()
            e1 = shape_data.observed_shape.e1
            e2 = shape_data.observed_shape.e2
            U[j,i] = np.sqrt(e1**2+e2**2)
            C[j,i] = shape_data.moments_sigma
            if C[j,i]<0:
                print 'negative size detected :/'
            else: 
                pass
            angles[j,i] = 0.5*np.arctan2(e2,e1)*180/np.pi
            if shape_data.error_message=='':
                pass
            else:
                print shape_data.error_message
                continue
    if ccd == 0:
        pass
    else:
        old_X = np.append(X,old_X)
        old_Y = np.append(Y,old_Y)
        old_U = np.append(U,old_U)
        old_C = np.append(C,old_C)
        old_angles = np.append(angles, old_angles)
        #        print old_X, old_Y
            # np.append
std2fwhm = 2.*np.sqrt(2.*np.log(2.))
#ind = np.lexsort((pixOx[:,0],pixOx[:,1]))
#ind = ind.astype(int)
#from scipy.stats import mode
#i,j = mode(pixOx,axis=0)
#j = j[0]
#i = j[0]
#j = j[1]
#f, axs = plt.subplots(2, 2, sharex='col', sharey='row',figsize=(10,12))
print np.shape(old_U)
old_V = np.zeros(np.shape(old_U))
print np.shape(old_V)
#for i,x in enumerate(np.linspace(0, width, nx)):
#    for j,y in enumerate(np.linspace(0, height, ny)):
#        image = psf.computeKernelImage(lsst.afw.geom.Point2D(x, y))
#        array = image.getArray()
#        galsim_image = galsim.Image(array)
#        shape_data = galsim.hsm.FindAdaptiveMom(galsim_image)
#        X[ccd,j,i] = x
#        Y[ccd,j,i] = y
#        e1 = shape_data.observed_shape.e1
#        e2 = shape_data.observed_shape.e2
#        U[ccd,j,i] = np.sqrt(e1**2+e2**2)
#        C[ccd,j,i] = shape_data.moments_sigma
#        angles[ccd,j,i] = 0.5*np.arctan2(e2,e1)*180/np.pi
#f, axs = plt.subplots(14,8, sharex='col', sharey='row',figsize=(40,40))
#for uh in range(8):
#    for hm in range(14):
#        axs[uh][hm].quiver(X[uh][hm], Y[uh][hm],U[uh][hm],V[uh][hm],angles=angles[uh][hm],headlength=0,headaxislength=0)
pixel_scale = r2arcs(pixel_scale)
plt.figure()
Q = plt.quiver(pixel_scale*(old_X-np.median(old_X)),pixel_scale*(old_Y-np.median(old_Y)),old_U,old_V,std2fwhm*old_C,angles=old_angles,headlength=0,headaxislength=0,scale=20,cmap='viridis')
qk = plt.quiverkey(Q, 0.9, 0.9, 0.05, r'$|e|=0.05$', labelpos='E',
                   coordinates='figure')
plt.title("whisker plot")
plt.colorbar(cmap='tab20b',label='FWHM')
#plt.axis('equal')

plt.xlabel('arcsec')
plt.ylabel('arcsec')
plt.savefig("/global/homes/h/husni/whisker.png")
plt.close()

#plt.figure(2)
#plt.plot(positions,'bo')
#plt.savefig("/global/homes/h/husni/positions.png")
#plt.close()
