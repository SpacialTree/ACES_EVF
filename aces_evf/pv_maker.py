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
import spectral_cube
from spectral_cube import SpectralCube

save_path = '/orange/adamginsburg/ACES/broadline_sources/EVFs/images/'

def make_position_list(amin, amax, step, unit):
    return (np.arange(amin.to(unit).value, amax.to(unit).value+1, step.value)*unit).to(u.deg)

def make_subcube():
    print('ok')

def main():
    basepath = '/orange/adamginsburg/ACES/mosaics/cubes/'
    cube_fn = f'{basepath}/CS_CubeMosaic.fits'
    mol = cube_fn.split('/')[-1].split('_')[0]
    cube = SpectralCube.read(cube_fn)
    
    # Latitude
    b_min = -0.27*u.deg
    b_max = 0.22*u.deg
    list_b = make_position_list(b_min, b_max, 0.5*u.arcmin, u.arcmin)

    # Longitude
    l_min = -0.59*u.deg
    l_max = 0.88*u.deg
    list_l = make_position_list(l_min, l_max, 0.5*u.arcmin, u.arcmin)

    for b in list_b:
        reg = regions.RectangleSkyRegion(center=SkyCoord((l_min+l_max)/2., b, frame='galactic'), width=1.5*u.deg, height=1*u.arcmin)
        #reg.to_pixel(WCS(head)).plot(ax=ax, edgecolor='red', facecolor='none')
        subcube = cube.subcube_from_regions([reg])
        pv_mean = subcube.mean(axis=1)
        pv_mean.save(f'{save_path}/{mol}_pv_b{round(b.value, 3)}_mean.fits')
        pv_max = subcube.max(axis=1)
        pv_max.save(f'{save_path}/{mol}_pv_b{round(b.value, 3)}_max.fits')

    for l in list_l:
        reg = regions.RectangleSkyRegion(center=SkyCoord(l, (b_min+b_max)/2., frame='galactic'), width=1*u.arcmin, height=0.5*u.deg)
        #reg.to_pixel(WCS(head)).plot(ax=ax, edgecolor='red', facecolor='none')
        subcube = cube.subcube_from_regions([reg])
        pv_mean = subcube.mean(axis=2)
        pv_mean.save(f'{save_path}/{mol}_pv_l{round(l.value, 3)}_mean.fits')
        pv_max = subcube.max(axis=2)
        pv_max.save(f'{save_path}/{mol}_pv_l{round(l.value, 3)}_max.fits')


if __name__ == "__main__":
    main()#sys.argv[1])