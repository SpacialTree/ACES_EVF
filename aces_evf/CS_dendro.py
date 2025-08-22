import numpy as np                          # For math and working with arrays
import matplotlib                           # For plotting! Lots of documentation online
import matplotlib.pyplot as plt
import astropy.io.fits as fits              # For importing the FITS image 
from astropy.table import Table, join       # Used for creating tables via astropy
from astropy.wcs import WCS                 # For adding coordinates to images
from astropy import units as u              # Helps with unit conversion!
from spectral_cube import SpectralCube      # Useful package for working with 3D data cubes
from astropy.coordinates import SkyCoord
from astropy.io import fits, ascii
import astrodendro
from astrodendro import dendrogram 
from astrodendro.analysis import PPStatistic
from astrodendro.structure import Structure
from matplotlib.patches import Rectangle
from regions import Region, RectanglePixelRegion, RectangleSkyRegion, Regions
from astropy.visualization import LogStretch, SqrtStretch, PercentileInterval, ImageNormalize
import glob
import os

cubepath = '/orange/adamginsburg/ACES/mosaics/cubes/'
savepath = '/orange/adamginsburg/ACES/broadline_sources/EVFs/dendro/'
file = 'CS21_CubeMosaic.fits'

# Grab the CS data cube:
cube = SpectralCube.read(cubepath+file)
# Calculate the beam area in deg^2 (gaussian area)
beam_omega_deg2 = (cube.header['BMAJ']*cube.header['BMIN']*np.pi)/(4*np.log(2))

# Convert beam parameters to get # of pixels in beam:
pix_area_deg2 = abs(cube.header['CDELT1']) * abs(cube.header['CDELT2'])
mpix = beam_omega_deg2/pix_area_deg2

# Dendro parameters --> currently set pretty high just to see if code runs
cube_noise = 0.1  #insert mad_std value of cube here
min_val = 100  # Minimum value factor 
min_sig = 100  # Minimum significance factor
scale = cube_noise   # CS Cube mad_std rms noise value

# Generate the dendrogram:
dend = astrodendro.Dendrogram.compute(cube.unmasked_data[:].value, 
                                      wcs = cube.wcs,
                                      min_value= min_val*scale, 
                                      min_delta= min_sig*scale, 
                                      min_npix=int(mpix))


# Save the dendrogram if you only want to generate it one time:
dend.save_to(savepath+'CS_HVCC_dendrogram_mv_'+str(min_val)+'_ms_'+str(min_sig)+'_mp_'+str(mpix)+'.fits')