from spectral_cube import SpectralCube
import regions
from regions import Regions
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import simple_norm
from astropy.coordinates import SkyCoord
import astropy.units as u

fn_CS = '/orange/adamginsburg/ACES/mosaics/cubes/CS21_CubeMosaic.fits'
cube_CS = SpectralCube.read(fn_CS)

regions = Regions.read('/blue/adamginsburg/savannahgramze/ACES_EVF/aces_evf/region_list.reg')

for reg in regions:
    subcube = cube_CS.subcube_from_regions([reg])
    subcube = subcube.with_spectral_unit(u.km/u.s, velocity_convention='radio')
    l = reg.center.galactic.l.value
    b = reg.center.galactic.b.value
    if l > 90:
        l = l - 360
    l = round(l, 2)

    subcube.write(f'/orange/adamginsburg/ACES/broadline_sources/EVFs/cubes/CS21_l{l}_b{b}.fits', overwrite=True)