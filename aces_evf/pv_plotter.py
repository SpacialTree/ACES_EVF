import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from pv_maker import make_position_list

# Path to where to save the PV diagrams
save_path = '/orange/adamginsburg/ACES/broadline_sources/EVFs/images/'

# Latitude extrema
B_MIN = -0.27 * u.deg
B_MAX = 0.22 * u.deg

# Longitude extrema
L_MIN = -0.59 * u.deg
L_MAX = 0.88 * u.deg

def plot_pv(pv, mol, pos, longitude=True, mean='mean', close=True):
    """
    Plot the PV diagram.

    Parameters
    ----------
    pv : SpectralCube
        SpectralCube object of the PV diagram.
    mol : str
        Molecule name.
    pos : Quantity
        Position of the PV diagram.
    longitude : bool
        If True, the PV diagram is along a longitude (l=pos). If False, the PV diagram is along a latitude (b=pos).
    """
    ww = WCS(pv.header)
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111, projection=ww)
    vmin = np.nanpercentile(pv.data, 1)
    vmax = np.nanpercentile(pv.data, 99)
    ax.imshow(pv.data, aspect='auto', vmin=vmin, vmax=vmax, cmap='inferno')

    ax11 = ax.coords[1]
    ax11.set_format_unit(u.km / u.s)
    ax.set_ylabel('Velocity (km/s)')

    if longitude:
        ax.set_xlabel('Latitude (deg)')
        ax.set_title(f'{mol} {mean} PV diagram at l={round(pos.value, 3)}')
        plt.tight_layout()
        plt.savefig(f'{save_path}/{mol}_{mean}_pv_l{round(pos.value, 3)}.png', overwrite=True)
    else:
        ax.set_xlabel('Longitude (deg)')
        ax.set_title(f'{mol} {mean} PV diagram at b={round(pos.value, 3)}')
        plt.tight_layout()
        plt.savefig(f'{save_path}/{mol}_{mean}_pv_b{round(pos.value, 3)}.png', overwrite=True)
    if close:
        plt.close()

def plot_all_pvb(opt='max'):
    list_b = make_position_list(B_MIN, B_MAX)
    for pos in list_b:
        pv_fn = f'{save_path}/CS21_pv_b{round(pos.value, 3)}_{opt}.fits'
        try:
            pv1 = fits.open(pv_fn)
        except FileNotFoundError:
            return False
        plot_pv(pv1[0], 'CS21', pos, longitude=False, mean=opt, close=True)

def plot_all_pvl(opt='max'):
    list_l = make_position_list(L_MIN, L_MAX)
    for pos in list_l:
        pv_fn = f'{save_path}/CS21_pv_l{round(pos.value, 3)}_{opt}.fits'
        try:
            pv1 = fits.open(pv_fn)
        except FileNotFoundError:
            return False
        plot_pv(pv1[0], 'CS21', pos, longitude=True, mean=opt, close=True)

def main():
    """
    Main function to make plots for each PV diagram.
    """

    plot_all_pvb(opt='max')
    plot_all_pvl(opt='max')
    plot_all_pvb(opt='mean')
    plot_all_pvl(opt='mean')

if __name__ == '__main__':
    main()
