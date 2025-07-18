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
from astropy import log
import os
import timeit
from astropy.table import Table
print('Running EVF Cutouts', flush=True)

#cube_mosaic_dir = '/orange/adamginsburg/ACES/mosaics/cubes/'
#save_dir = '/orange/adamginsburg/ACES/broadline_sources/EVFs/cubes/'
# /home/savannahgramze/blue_adam/savannahgramze/ACES_EVF/aces_evf
#regs = Regions.read('/blue/adamginsburg/savannahgramze/ACES_EVF/aces_evf/EVF_reg_list.reg')
#line_list = ['CS21']#['H13CN', 'H13COp', 'SiO21', 'SO32', 'SO21', 'HN13C', 'HC3N', 'CS21']
#H13CN, H13CO+, SiO, HCO+, HNCO, SO(3-2), SO(2-1), HN13C, HC3N, CS(2-1)
# test

def main(line_list=['CS21'], mode=False): 
    save_dir = '/orange/adamginsburg/ACES/broadline_sources/EVFs/cubes/'
    for line in line_list:
        line_dir = os.path.join(save_dir, line)
        if not os.path.exists(line_dir):
            print(f'Creating directory {line_dir}', flush=True)
            os.makedirs(line_dir)

    #regs = Regions.read('/blue/adamginsburg/savannahgramze/ACES_EVF/aces_evf/EVF_reg_list.reg')
    tbl = Table.read('/blue/adamginsburg/savannahgramze/ACES_EVF/aces_evf/Filtered_EVFs_table.ecsv')
    if not mode: 
        print('Processing all regions individually', flush=True)
        if len(line_list) > 1:
            for row in tbl: 
                l = row['l']
                b = row['b']
                reg = regions.RectangleSkyRegion(SkyCoord(l, b, unit='deg', frame='galactic'), width=row['deltal']*u.deg, height=row['deltab']*u.deg)
                l = round(l, 3)
                b = round(b, 3)

                for line in line_list:
                    line_dir = os.path.join(save_dir, line)
                    fn = os.path.join(cube_mosaic_dir, f'{line}_CubeMosaic.fits')
                    cube = SpectralCube.read(fn, use_dask=False)
                    

                    # Create a subcube from the region
                    start = timeit.default_timer()
                    subcube = cube.subcube_from_regions([reg])
                    subcube = subcube.with_spectral_unit(u.km/u.s, velocity_convention='radio')
                    end = timeit.default_timer()
                    print(f'Subcube creation took {end - start} seconds', flush=True)
                    
                    write_fn = os.path.join(line_dir, f'EVF_{row["ID Number"]}_{line}_l{l}_b{b}.fits')
                    print(f'Writing {line} l{l} b{b} to {write_fn}', flush=True)
                    start = timeit.default_timer()
                    subcube.write(write_fn, overwrite=True)
                    end = timeit.default_timer()
                    print(f'Writing took {end - start} seconds', flush=True)
        elif len(line_list) == 1:
            for line in line_list:
                line_dir = os.path.join(save_dir, line)
                fn = os.path.join(cube_mosaic_dir, f'{line}_CubeMosaic.fits')
                print(f"Reading cube: {fn}", flush=True)
                cube = SpectralCube.read(fn, use_dask=False)
                #tbl = Table.read('/blue/adamginsburg/savannahgramze/ACES_EVF/aces_evf/Filtered_EVFs_table.ecsv')

                print('Beginning to loop through individual fields to make and write subcubes.', flush=True)
                for row in tbl: 
                    l = row['l']#reg.center.galactic.l.value
                    b = row['b']#reg.center.galactic.b.value
                    reg = regions.RectangleSkyRegion(SkyCoord(l, b, unit='deg', frame='galactic'), width=row['deltal']*u.deg, height=row['deltab']*u.deg)
                    l = round(l, 3)
                    b = round(b, 3)

                    # Create a subcube from the region
                    start = timeit.default_timer()
                    subcube = cube.subcube_from_regions([reg])
                    subcube = subcube.with_spectral_unit(u.km/u.s, velocity_convention='radio')
                    end = timeit.default_timer()
                    print(f'Subcube creation took {end - start} seconds', flush=True)

                    write_fn = os.path.join(line_dir, f'EVF_{row["ID Number"]}_{line}_l{l}_b{b}.fits')
                    print(f'Writing {line} l{l} b{b} to {write_fn}', flush=True)
                    start = timeit.default_timer()
                    subcube.write(write_fn, overwrite=True)
                    end = timeit.default_timer()
                    print(f'Writing took {end - start} seconds', flush=True)
        else:
            print('No lines to process', flush=True)
            exit(1)
    else: 
        print('Processing all regions together', flush=True)
        for line in line_list:
            line_dir = os.path.join(save_dir, line)
            fn = os.path.join(cube_mosaic_dir, f'{line}_CubeMosaic.fits')
            cube = SpectralCube.read(fn, use_dask=False)

            # Create a subcube from the region list
            print(f'Creating subcube for {line} from all regions', flush=True)
            start = timeit.default_timer()
            subcube = cube.subcube_from_regions(regs)
            subcube = subcube.with_spectral_unit(u.km/u.s, velocity_convention='radio')
            end = timeit.default_timer()
            print(f'Subcube creation took {end - start} seconds', flush=True)

            print(f'Writing {line} to {os.path.join(line_dir, f"{line}_EVF_masked.fits")}', flush=True)
            start = timeit.default_timer()
            subcube.write(os.path.join(line_dir, f'{line}_EVF_masked.fits'), overwrite=True)
            end = timeit.default_timer()
            print(f'Writing took {end - start} seconds', flush=True)

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-l", "--lines", dest="lines", default='CS21',
                      help="comma-separated list of lines to process", metavar="lines")
    parser.add_option('-m', '--big', dest='big', action='store_true',# default=True,
                      help='On=all regions together, Off=separate regions')
    (options, args) = parser.parse_args()
    lines = options.lines.split(',')
    mode = bool(options.big)
    print(f'Processing lines: {lines} with mode: {mode}', flush=True)

    cube_mosaic_dir = '/orange/adamginsburg/ACES/mosaics/cubes/'
    for line in lines:
        if not os.path.exists(os.path.join(cube_mosaic_dir, f'{line}_CubeMosaic.fits')):
            print(f'Cube for {line} does not exist in {cube_mosaic_dir}. Please check the line name.', flush=True)
            exit(1)

    main(lines, mode)

#for line in line_list:
#    line_dir = os.path.join(save_dir, line)
#    if not os.path.exists(line_dir):
#        os.makedirs(line_dir)
#
#for reg in regions:
#    l = reg.center.galactic.l.value
#    b = reg.center.galactic.b.value
#    if l > 90:
#        l = l - 360
#    l = round(l, 3)
#    b = round(b, 3)
#
#    for line in line_list:
#        line_dir = os.path.join(save_dir, line)
#        fn = os.path.join(cube_mosaic_dir, f'{line}_CubeMosaic.fits')
#        cube = SpectralCube.read(fn)
#
#        # Create a subcube from the region
#        subcube = cube.subcube_from_regions([reg])
#        subcube = subcube.with_spectral_unit(u.km/u.s, velocity_convention='radio')
#        subcube.write(os.path.join(line_dir, f'{line}_l{l}_b{b}.fits'), overwrite=True)

#large_lines = ['HNCO_7m12mTP', 'HCOP_noTP']
#
#for line in large_lines:
#    line_dir = os.path.join(save_dir, line)
#    if not os.path.exists(line_dir):
#        os.makedirs(line_dir)
#
#for reg in regions:
#    l = reg.center.galactic.l.value
#    b = reg.center.galactic.b.value
#    if l > 90:
#        l = l - 360
#    l = round(l, 2)
#    b = round(b, 2)
#
#    for line in large_lines:
#        line_dir = os.path.join(save_dir, line)
#        fn = os.path.join(cube_mosaic_dir, f'{line}_CubeMosaic.fits')
#        cube = SpectralCube.read(fn)
#
#        # Create a subcube from the region
#        subcube = cube.subcube_from_regions([reg])
#        subcube = subcube.with_spectral_unit(u.km/u.s, velocity_convention='radio')
#        subcube.write(os.path.join(line_dir, f'{line}_l{l}_b{b}.fits'), overwrite=True)
#

