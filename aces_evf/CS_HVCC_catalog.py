# This is a script for generating a catalog from the CS dendrogram:

import numpy as np                          # For math and working with arrays
import astropy.io.fits as fits              # For importing the FITS image 
from astropy.table import Table, join       # Used for creating tables via astropy
from astropy.wcs import WCS                 # For adding coordinates to images
from astropy import units as u              # Helps with unit conversion!
from astropy.coordinates import SkyCoord
from astropy.io import fits, ascii
import astrodendro
from astrodendro import dendrogram 
from astrodendro.analysis import PPStatistic
from astrodendro.structure import Structure


# Set up file paths: 
tabfile = # Add by-eye ID table path+file name here
cubefile = # Add high resolution CS cube path+file name here
dendfile = # Add dendrogram path+file name here
catfile = # Save Path + filename of the catalog 

# Grab the CS data cube:
cube = SpectralCube.read(cubefile)
wcs = cube.wcs

# Calculate the beam area in deg^2 (gaussian area)
beam_omega_deg2 = (cube.header['BMAJ']*cube.header['BMIN']*np.pi)/(4*np.log(2))

# Convert beam parameters to get # of pixels in beam:
pix_area_deg2 = abs(cube.header['CDELT1']) * abs(cube.header['CDELT2'])
mpix = beam_omega_deg2/pix_area_deg2

# Read in dendrogram file from high res CS cube
dend = astrodendro.Dendrogram.load_from(dendfile)

# Create the catalog from the dendrogram:
Metadata = {}

# Define the pixel size in arcseconds:
PixelAreaArcsec = 3600. * abs(cube.header['CDELT1']) * 3600. * abs(cube.header['CDELT2'])
Metadata['wcs'] = wcs
# Define the pixel units in the image input to the dendrogram:
Metadata['data_unit'] = u.Jy / u.beam
# Define the pixel size (on one side):
Metadata['spatial_scale'] =  PixelAreaArcsec**0.5 * u.arcsec
Metadata['velocity_scale'] = cube.header['CDELT3'] * u.m/u.s
# Define the beam parameters (since pixel units are Jy/beam):
Metadata['beam_major'] = (3600 * cube.header['BMAJ'])*u.arcsec
Metadata['beam_minor'] = (3600 * cube.header['BMIN'])*u.arcsec

catalog = astrodendro.ppv_catalog(dend, Metadata, verbose=False)

catalog['idx'] = [structure.idx for structure in dend]
indices = [dend[index].indices() for i, index in enumerate(catalog['idx'])]

catalog['vmax_pix'], catalog['vmin_pix'] = zip(*[(max(idx[0]), min(idx[0])) for idx in indices])
catalog['lmax_pix'], catalog['lmin_pix'] = zip(*[(max(idx[2]), min(idx[2])) for idx in indices])
catalog['bmax_pix'], catalog['bmin_pix'] = zip(*[(max(idx[1]), min(idx[1])) for idx in indices])

#convert to WCS
sky_coords, velocity_coords = wcs.pixel_to_world(catalog['lmax_pix'], catalog['bmax_pix'], catalog['vmax_pix'])
catalog['lmax_u'], catalog['bmax_u'] = sky_coords.l.deg, sky_coords.b.deg
catalog['vmax_u'] = velocity_coords
sky_coords, velocity_coords = wcs.pixel_to_world(catalog['lmin_pix'], catalog['bmin_pix'], catalog['vmin_pix'])
catalog['lmin_u'], catalog['bmin_u'] = sky_coords.l.deg, sky_coords.b.deg
catalog['vmin_u'] = velocity_coords

catalog['lmax_u'] = np.where(catalog['lmax_u'] > 180, catalog['lmax_u'] - 360, catalog['lmax_u']) 
catalog['lmin_u'] = np.where(catalog['lmin_u'] > 180, catalog['lmin_u'] - 360, catalog['lmin_u']) 


catalog['del_v']= catalog['vmax_u']-catalog['vmin_u']
catalog['del_l']= catalog['lmax_u']-catalog['lmin_u']
catalog['del_b']= catalog['bmax_u']-catalog['bmin_u']
catalog=catalog[(catalog['del_l']!= 0)] 

catalog['verticality'] = np.divide(np.abs(catalog['del_v']), np.abs(catalog['del_l'])) 
catalog['verticality'].unit = u.m/(u.s*u.deg)

peakl,peakb,peakv = [dend[index].get_peak()[0][2] for i, index in enumerate(catalog['idx'])], [dend[index].get_peak()[0][1] for i,index in enumerate(catalog['idx'])], [dend[index].get_peak()[0][0] for i,index in enumerate(catalog['idx'])]
peak_sky_coords, peak_velocity_coords = wcs.pixel_to_world(peakl, peakb, peakv)
catalog['l_peak'], catalog['b_peak'] = peak_sky_coords.l.deg, peak_sky_coords.b.deg
catalog['v_peak']= peak_velocity_coords

# Create a mask that removes sources based on verticality, l extent and v extent:
cat_mask = (((catalog['lmax_pix']-catalog['lmin_pix'])>=2.) 
            & (catalog['verticality'] >= 1000000)
            & ((catalog['lmax_pix']-catalog['lmin_pix'])<50.)
            & (np.abs(catalog['del_v'] > 10000)))

catalog = catalog[cat_mask]

catalog.write(catfile, overwrite=True)