import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import lsst.daf.persistence
import galsim

# Create a data butler which provides access to a data repository.
butler = lsst.daf.persistence.Butler("/project/projectdirs/lsst/zuntz/reductions/run2")

# I don't know if this is the proper way to get hold of a list of which visits (individual pointing) and CCDs (chips in the focal plane) 
# are in the data set, but it is the think I could find.
mapper=butler.repository.mappers()[0]
expmapping=mapper.exposures['calexp']
#Get the list of pairs of (visit,ccd) in the database
ccd_exposures = expmapping.lookup(("visit", "ccd"), {})
#Cut down to the ones that actually exist.  I think the other ones must not have processed correctly or something like that
ccd_exposures = [(visit,ccd) for (visit,ccd) in ccd_exposures if butler.datasetExists("calexp", visit=visit, ccd=ccd)]


#For this example we just look at the first chip
visit, ccd = ccd_exposures[0]

#Ask the Butler for the calibration exposure  for this visit.
calexp = butler.get("calexp", visit=visit, ccd=ccd, immediate=True)

# Get the  width and height of the image and an object representing the PSF variation across the image.
width = calexp.getWidth()
height = calexp.getHeight()
psf = calexp.getPsf()
print width/height
# As an example let's make a 200 x 100 image of the variation of the PSF ellipticity 
# across the chip
nx = 100
ny = int(nx * (height/width))

#Maps we will fill in with the results
psf_g1_map = np.zeros((ny,nx))
psf_g2_map = np.zeros((ny,nx))

#Loop through pixels in the chip.  We just take nx and ny samples across the image
# not every single pixel
for i,x in enumerate(np.linspace(0, width, nx)):
    for j,y in enumerate(np.linspace(0, height, ny)):
        # Compute an LSST-format image of the PSF at this point
        image = psf.computeKernelImage(lsst.afw.geom.Point2D(x, y))
        # Extract the numpy array from it
        array = image.getArray()

        #Make this into a galsim image.
        galsim_image = galsim.Image(array)

        # Measure the shape of the PSF image and fill in our map pixel.
        shape_data = galsim.hsm.FindAdaptiveMom(galsim_image,strict=False)
        psf_g1_map[j,i] =shape_data.moments_sigma

# Save a plot of the image
plt.figure()
plt.imshow(psf_g1_map, interpolation='nearest')
plt.title("g1")
plt.savefig("/global/homes/h/husni/sigma_map.png")
plt.close()
