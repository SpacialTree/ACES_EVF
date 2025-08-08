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

# Function to help with cross-matching dendrogram & by-eye IDs
def is_coordinate_in_3d_range(point, min_corner, max_corner):
    """
    Checks if a 3D coordinate is within a specified 3D range.

    Args:
        point (tuple): A tuple (x, y, z) representing the 3D coordinate.
        min_corner (tuple): A tuple (min_x, min_y, min_z) representing the minimum corner of the range.
        max_corner (tuple): A tuple (max_x, max_y, max_z) representing the maximum corner of the range.

    Returns:
        bool: True if the point is within the range, False otherwise.
    """
    x, y, z = point
    min_x, min_y, min_z = min_corner
    max_x, max_y, max_z = max_corner

    return (min_x <= x <= max_x and
            min_y <= y <= max_y and
            min_z <= z <= max_z)

cubefile = # high res CS cube path+file
tabfile = # by-eye table path+file
catfile = # catalog path+file
match_cat = # Save path+file for catalog corresponding to only the by-eye EVFs

# Grab the CS data cube:
cube = SpectralCube.read(cubefile)  
wcs = cube.wcs

# Grab the EVF by-eye ID table:
tbl = Table(ascii.read(tabfile, format='ecsv'))

# Get EVF coords in pixels:
center_sky = SkyCoord(tbl['l'], tbl['b'], frame = 'galactic')

tbl['vert'] = np.divide(np.abs(tbl['deltaV']), np.abs(tbl['deltal'])) 
tbl['x_pix'], tbl['y_pix'] = center_sky.to_pixel(wcs)
tbl['v_pix'] =(((tbl['V_LSR'].data- cube.header['CRVAL3'])/(cube.header['CDELT3'] )) + cube.header['CRPIX3'] )#.astype(int)

# Currently rounding values and not using int:
tbl['x_pix'] = [round(item) for item in tbl['x_pix']]
tbl['y_pix'] = [round(item) for item in tbl['y_pix']]
tbl['v_pix'] = [round(item) for item in tbl['v_pix']]


#Read in the dendro catalog:
cat = Table.read(catfile)

# Cross match dendrogram catalog with by-eye detections:
d = []
for c, coord in enumerate(tbl['l']):
    pix_coords = (tbl['x_pix'][c], tbl['y_pix'][c], tbl['v_pix'][c])

    x=[]
    area= []
    # Check to see if by-eye coords are in dendrogram extents:
    for s, struc in enumerate(cat['idx']):
        cond = is_coordinate_in_3d_range(pix_coords,
                                        (cat['lmin_pix'][s], cat['bmin_pix'][s],  cat['vmin_pix'][s]), 
                                         (cat['lmax_pix'][s], cat['bmax_pix'][s], cat['vmax_pix'][s])) 
        # Only keep structures that match with by-eye detections:
        if cond == True:
            x.append(struc)
            area.append(cat['area_exact'][s])
        
    if x == []:
        d.append(np.nan)
    else:
        x = np.array(x)
        area = np.array(area)
        max_x = x[area == np.nanmax(area)]
        
        d.extend(x) #Use append(max_x[0] instead if you just want the largest structure containing the given coordinate

d = np.array(d)
cat.add_index('idx')
cat = cat.loc[np.unique(d[~np.isnan(d)])]

cat.write(match_cat, overwrite=True)