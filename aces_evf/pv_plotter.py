
from pv_maker import make_position_list

# Path to where to save the PV diagrams
save_path = '/orange/adamginsburg/ACES/broadline_sources/EVFs/images/'

# Latitude extrema
B_MIN = -0.27 * u.deg
B_MAX = 0.22 * u.deg

# Longitude extrema
L_MIN = -0.59 * u.deg
L_MAX = 0.88 * u.deg

def plot_all_pvb(opt='max'):
    list_b = make_position_list(B_MIN, B_MAX)
    for pos in list_b:
        pv_fn = f'{save_path}/CS21_pv_b{round(pos.value, 3)}_{opt}.fits'
        pv1 = fits.open(pv_fn)
        plot_pv(pv1[0], 'CS21', pos, longitude=False, mean=opt, close=True)

def plot_all_pvl(opt='max'):
    list_l = make_position_list(L_MIN, L_MAX)
    for pos in list_l:
        pv_fn = f'{save_path}/CS21_pv_l{round(pos.value, 3)}_{opt}.fits'
        pv1 = fits.open(pv_fn)
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
