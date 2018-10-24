import sys, os, math, galsim
from astropy.io import fits
import numpy as np


def main(argv):
    gal_flux = 1.e5    # total counts on the image
    gal_sigma = 2.     # arcsec
    psf_sigma = 1.     # arcsec
    pixel_scale = 0.2  # arcsec / pixel
    noise, image, result_moments = [[0 for _ in range(100)] for _ in range(3)]
    noise_array = np.ones(100)*np.random.rand(100)*2
    for i in range(100):
        noise[i] = 30.*noise_array[i]

    # Define the galaxy profile
    gal = galsim.Gaussian(flux=gal_flux, sigma=gal_sigma)
    gal_im = gal.drawImage(scale=pixel_scale)

    # Define the PSF profile
    psf = galsim.Gaussian(flux=1., sigma=psf_sigma) # PSF flux should always = 1
    psf_im = psf.drawImage(scale=pixel_scale)

    # Final profile is the convolution of these
    final = galsim.Convolve([gal, psf])

    # Draw the image with a particular pixel scale, given in arcsec/pixel.
    for i in range(100):
        image[i] = final.drawImage(scale=pixel_scale)

    # Add Gaussian noise to the image with specified sigma
    for i in range(100):
        image[i].addNoise(galsim.GaussianNoise(sigma=noise[i]))

    # Write the image to a file
    for i in range(100):
        file_name = 'test'+str(i)+'.fits'
        image[i].write(file_name)
    for i in range(100):
        result_moments[i] = image[i].FindAdaptiveMom()
    result_shear = galsim.hsm.EstimateShear(gal_im,psf_im)

    return result_moments, result_shear
if __name__ == "__main__":
    main(sys.argv)
