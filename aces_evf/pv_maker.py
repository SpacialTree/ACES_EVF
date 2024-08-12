from astropy.io import fits
import numpy as np
import os
import sys
import glob
import pvextractor
from pvextractor import Path, extract_pv_slice
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from astropy.visualization import simple_norm

def main():
    cube_fn = '/orange/adamginsburg/ACES/mosaics/cubes/CS_CubeMosaic.fits'
    mol = cube_fn.split('/')[-1].split('_')[0]
    cube = fits.open(cube_fn)
    wcs = WCS(cube[0].header)
    cube_data = cube[0].data
    
    # Latitude
    b_min = -0.27*u.deg
    b_max = 0.22*u.deg
    list_b = np.arange(b_min.to(u.arcmin).value, b_max.to(u.arcmin).value+1, 1)*u.arcmin
    list_b = list_b.to(u.deg)

    # Longitude
    l_min = -0.59*u.deg
    l_max = 0.88*u.deg
    list_l = np.arange(l_min.to(u.arcmin).value, l_max.to(u.arcmin).value+1, 1)*u.arcmin
    list_l = list_l.to(u.deg)

    # Make PV diagrams across all longitudes
    for i in range(len(list_l)):
        l = list_l[i]
        c1 = SkyCoord(l, b_min, frame='galactic')
        c2 = SkyCoord(l, b_max, frame='galactic')
        path = Path([c1, c2], width=1*u.arcmin)
        pv = extract_pv_slice(cube_data, path, wcs)
        pv.writeto(f'pv_{mol}_l{l.value}.fits', overwrite=True)

    # Make PV diagrams across all latitudes
    for i in range(len(list_b)):
        b = list_b[i]
        c1 = SkyCoord(l_min, b, frame='galactic')
        c2 = SkyCoord(l_max, b, frame='galactic')
        path = Path([c1, c2], width=1*u.arcmin)
        pv = extract_pv_slice(cube_data, path, wcs)
        pv.writeto(f'pv_{mol}_b{b.value}.fits', overwrite=True)

    

if __name__ == "__main__":
    main()#sys.argv[1])