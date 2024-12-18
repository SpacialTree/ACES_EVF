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
import regions
from regions import Regions
import time

# Path to where to save the PV diagrams
save_path = '/orange/adamginsburg/ACES/broadline_sources/EVFs/images/HNCO/'

# Latitude extrema
B_MIN = -0.27 * u.deg
B_MAX = 0.22 * u.deg

# Longitude extrema
L_MIN = -0.59 * u.deg
L_MAX = 0.88 * u.deg


def make_position_list(amin, amax, step=0.5 * u.arcmin, unit=u.arcmin):
    """
    Make a list of Sky positions from amin to amax with a given step size and unit.
    Ideally used to make a list of latitudes or longitudes for the centers of region cutouts.

    Parameters
    ----------
    amin : Quantity
        Minimum value of the list.
    amax : Quantity
        Maximum value of the list.
    step : Quantity
        Step size between values.
    unit : astropy.unit
        Unit of the values.
    """

    return (np.arange(amin.to(unit).value, amax.to(unit).value + 1, step.value) * unit).to(u.deg)


def make_pv_b(cube, latitude, mol):
    """
    Make a PV diagram along a given latitude.

    Parameters
    ----------
    cube : SpectralCube
        SpectralCube object of the data.
    latitude : Quantity
        Latitude of the center of the region.
    mol : str
        Molecule name.
    """

    reg = regions.RectangleSkyRegion(center=SkyCoord((L_MIN + L_MAX) / 2., latitude, frame='galactic'), width=1.5 * u.deg, height=1 * u.arcmin)
    subcube = cube.subcube_from_regions([reg])
    subcube.allow_huge_operations = True

    t = time.process_time()
    subcube.rechunk(save_to_tmp_dir=True)
    print(f'Time to rechunk: {time.process_time() - t}')

    with subcube.use_dask_scheduler('threads', num_workers=4):

        t = time.process_time()
        pv_mean = subcube.mean(axis=1)
        print(f'Time to calculate mean: {time.process_time() - t}')

        t = time.process_time()
        pv_mean.write(f'{save_path}/{mol}_pv_b{round(latitude.value, 3)}_mean.fits', overwrite=True)
        print(f'Time to write mean: {time.process_time() - t}')

        t = time.process_time()
        pv_max = subcube.max(axis=1)
        print(f'Time to calculate max: {time.process_time() - t}')

        t = time.process_time()
        pv_max.write(f'{save_path}/{mol}_pv_b{round(latitude.value, 3)}_max.fits', overwrite=True)
        print(f'Time to write max: {time.process_time() - t}')

    print('Done with b =', latitude)


def make_pv_l(cube, longitude, mol):
    """
    Make a PV diagram along a given longitude.

    Parameters
    ----------
    cube : SpectralCube
        SpectralCube object of the data.
    longitude : Quantity
        Longitude of the center of the region.
    mol : str
        Molecule name.
    """

    reg = regions.RectangleSkyRegion(center=SkyCoord(longitude, (B_MIN + B_MAX) / 2., frame='galactic'), width=1 * u.arcmin, height=0.5 * u.deg)
    subcube = cube.subcube_from_regions([reg])
    subcube.allow_huge_operations = True

    t = time.process_time()
    subcube.rechunk(save_to_tmp_dir=True)
    print(f'Time to rechunk: {time.process_time() - t}')

    with subcube.use_dask_scheduler('threads', num_workers=4):

        t = time.process_time()
        pv_mean = subcube.mean(axis=2)
        print(f'Time to calculate mean: {time.process_time() - t}')

        t = time.process_time()
        pv_mean.write(f'{save_path}/{mol}_pv_l{round(longitude.value, 3)}_mean.fits', overwrite=True)
        print(f'Time to write mean: {time.process_time() - t}')

        t = time.process_time()
        pv_max = subcube.max(axis=2)
        print(f'Time to calculate max: {time.process_time() - t}')

        t = time.process_time()
        pv_max.write(f'{save_path}/{mol}_pv_l{round(longitude.value, 3)}_max.fits', overwrite=True)
        print(f'Time to write max: {time.process_time() - t}')

    print('Done with l =', longitude)

def make_pv(cube, position, mol, axis):
    """ 
    Make a PV diagram along a given axis.

    Parameters
    ----------
    cube : SpectralCube
        SpectralCube object of the data.
    position : Quantity
        Position of the center of the region.
    mol : str
        Molecule name.
    axis : int
        Axis along which to make the PV diagram. 
        1 for Galactic Latitude "b", 2 for Galactic Longitude "l".
    """

    if axis == 1:
        pv_axis = 'b'
        reg = regions.RectangleSkyRegion(center=SkyCoord((L_MIN + L_MAX) / 2., position, frame='galactic'), width=1.5 * u.deg, height=1 * u.arcmin)
    elif axis == 2:
        pv_axis = 'l'
        reg = regions.RectangleSkyRegion(center=SkyCoord(position, (B_MIN + B_MAX) / 2., frame='galactic'), width=1 * u.arcmin, height=0.5 * u.deg)
    else:
        print('Invalid axis')
        return

    subcube = cube.subcube_from_regions([reg])
    subcube.allow_huge_operations = True

    #t = time.process_time()
    #subcube.rechunk(save_to_tmp_dir=True)
    #print(f'Time to rechunk: {time.process_time() - t}')

    with subcube.use_dask_scheduler('threads', num_workers=4):
        
        t = time.process_time()
        pv_mean = subcube.mean(axis=axis)
        print(f'Time to calculate mean: {time.process_time() - t}')

        t = time.process_time()
        pv_mean.write(f'{save_path}/{mol}_pv_{pv_axis}{round(position.value, 3)}_mean.fits', overwrite=True)
        print(f'Time to write mean: {time.process_time() - t}')

        t = time.process_time()
        pv_max = subcube.max(axis=axis)
        print(f'Time to calculate max: {time.process_time() - t}')

        t = time.process_time()
        pv_max.write(f'{save_path}/{mol}_pv_{pv_axis}{round(position.value, 3)}_max.fits', overwrite=True)
        print(f'Time to write max: {time.process_time() - t}')

    print(f'Done with {pv_axis} =', position)

def make_pv_mol(cube_fn):
    """
    Make PV diagrams for a given molecule.

    Parameters
    ----------
    cube_fn : str
        Path to the fits file of the cube.
    """

    t = time.process_time()
    cube = SpectralCube.read(cube_fn, use_dask=True)
    print(f'Time to read cube: {time.process_time() - t}')
    #cube.use_dask_scheduler('threads', num_workers=4)
    mol = cube_fn.split('/')[-1].split('_')[0]

    list_b = make_position_list(B_MIN, B_MAX)
    list_l = make_position_list(L_MIN, L_MAX)

    print('Starting PV diagrams for', mol)
    for latitude in list_b:
        print('Making PV diagram at b =', latitude)
        make_pv(cube, latitude, mol, 1)

    for longitude in list_l:
        print('Making PV diagram at l =', longitude)
        make_pv(cube, longitude, mol, 2)


def main():
    basepath = '/orange/adamginsburg/ACES/mosaics/cubes/'
    #cube_fn_CS = f'{basepath}/CS21_CubeMosaic.fits'
    #make_pv_mol(cube_fn_CS)

    print('Starting HNCO')
    cube_fn_HNCO = f'{basepath}/HNCO_7m12mTP_CubeMosaic.fits'
    make_pv_mol(cube_fn_HNCO)


if __name__ == "__main__":
    main()#sys.argv[1])