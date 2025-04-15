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
import os

cube_mosaic_dir = '/orange/adamginsburg/ACES/mosaics/cubes/'
save_dir = '/orange/adamginsburg/ACES/broadline_sources/EVFs/cubes/'
# /home/savannahgramze/blue_adam/savannahgramze/ACES_EVF/aces_evf
regions = Regions.read('/blue/adamginsburg/savannahgramze/ACES_EVF/aces_evf/EVF_reg_list.reg')
line_list = ['H13CN', 'H13COp', 'SiO21', 'SO32', 'SO21', 'HN13C', 'HC3N', 'CS21']
#H13CN, H13CO+, SiO, HCO+, HNCO, SO(3-2), SO(2-1), HN13C, HC3N, CS(2-1)

for line in line_list:
    line_dir = os.path.join(save_dir, line)
    if not os.path.exists(line_dir):
        os.makedirs(line_dir)

for reg in regions:
    l = reg.center.galactic.l.value
    b = reg.center.galactic.b.value
    if l > 90:
        l = l - 360
    l = round(l, 2)
    b = round(b, 2)

    for line in line_list:
        line_dir = os.path.join(save_dir, line)
        fn = os.path.join(cube_mosaic_dir, f'{line}_CubeMosaic.fits')
        cube = SpectralCube.read(fn)

        # Create a subcube from the region
        subcube = cube.subcube_from_regions([reg])
        subcube = subcube.with_spectral_unit(u.km/u.s, velocity_convention='radio')
        subcube.write(os.path.join(line_dir, f'{line}_l{l}_b{b}.fits'), overwrite=True)

large_lines = ['HNCO_7m12mTP', 'HCOP_noTP']

for line in large_lines:
    line_dir = os.path.join(save_dir, line)
    if not os.path.exists(line_dir):
        os.makedirs(line_dir)

for reg in regions:
    l = reg.center.galactic.l.value
    b = reg.center.galactic.b.value
    if l > 90:
        l = l - 360
    l = round(l, 2)
    b = round(b, 2)

    for line in large_lines:
        line_dir = os.path.join(save_dir, line)
        fn = os.path.join(cube_mosaic_dir, f'{line}_CubeMosaic.fits')
        cube = SpectralCube.read(fn)

        # Check if the cube is empty
        if cube.size == 0:
            print(f"Cube {fn} is empty. Skipping...")
            continue

        # Check if the region is empty
        if reg.size == 0:
            print(f"Region {reg} is empty. Skipping...")
            continue

        # Create a subcube from the region
        subcube = cube.subcube_from_regions([reg])
        subcube = subcube.with_spectral_unit(u.km/u.s, velocity_convention='radio')
        subcube.write(os.path.join(line_dir, f'{line}_l{l}_b{b}.fits'), overwrite=True)
