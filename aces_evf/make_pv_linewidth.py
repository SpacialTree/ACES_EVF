#import libraries
import os, glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.colors as cm
from astropy.io import fits, ascii
from astropy.wcs import WCS
from astropy import units as u
from matplotlib.gridspec import GridSpec
import astropy.io.fits as pyfits
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel
from astropy.convolution import convolve_fft, Gaussian2DKernel
import matplotlib.colors as mc
from astropy.table import Table, vstack
from matplotlib.patches import FancyArrowPatch
from matplotlib.backends.backend_pdf import PdfPages
from spectral_cube import SpectralCube
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

#silencing slowdown warnings
import warnings
from spectral_cube.utils import SpectralCubeWarning
warnings.filterwarnings(action='ignore', category=SpectralCubeWarning,
                        append=True)


#write paths
basepath = '/orange/adamginsburg/ACES/broadline_sources/EVFs/'
CS_cubes_path   = basepath + 'cubes/CS21/'
SiO_cubes_path  = basepath + 'cubes/SiO21/'
HC3N_cubes_path = basepath + 'cubes/HC3N/'


##Open the filted EVF table 
EVF_tab = Table.read('Filtered_EVFs_table.ecsv')
EVF_v_cent  = EVF_tab['V_LSR']
EVF_v_width = EVF_tab['max_v'] - EVF_tab['min_v']


def make_pv_plots(cube,ID_num,EVF_filename):
    LV = cube.mean(axis=1)
    BV = cube.mean(axis=2)

    pyfits.writeto(basepath + 'pvdiagrams/EVF_{}_{}_LV.fits'.format(ID_num, EVF_filename), LV.value, LV.header, overwrite=True)
    pyfits.writeto(basepath + 'pvdiagrams/EVF_{}_{}_BV.fits'.format(ID_num, EVF_filename), BV.value, BV.header, overwrite=True)


def single_gauss(x,mu,sig, A):
    return A*np.exp(-(x-mu)**2/2/sig**2)

def double_gauss(x,mu1,mu2,A1,A2,sig1,sig2):
    return single_gauss(x,mu1,sig1,A1)+single_gauss(x,mu2,sig2,A2)

def make_pv_plots(cube,ID_num,EVF_filename):
    LV = cube.mean(axis=1)
    BV = cube.mean(axis=2)

    pyfits.writeto('./pvdiagrams/EVF_{}_{}_LV.fits'.format(ID_num, EVF_filename), LV.value, LV.header, overwrite=True)
    pyfits.writeto('./pvdiagrams/EVF_{}_{}_BV.fits'.format(ID_num, EVF_filename), BV.value, BV.header, overwrite=True)


def single_gauss(x,mu,sig, A):
    return A*np.exp(-(x-mu)**2/2/sig**2)

def double_gauss(x,mu1,mu2,A1,A2,sig1,sig2):
    return single_gauss(x,mu1,sig1,A1)+single_gauss(x,mu2,sig2,A2)


def fit_gaussians(cube, ID_num,EVF_filename):
    tab_ind = np.where(EVF_tab['ID Number'] == int(ID_num))[0][0]
    subcube = cube.spectral_slab(EVF_tab['min_v'][tab_ind]* u.km / u.s, EVF_tab['max_v'][tab_ind]* u.km / u.s)
    avg_spectrum = cube.mean(axis=(1,2))
    avg_subcube_spectrum = subcube.mean(axis=(1,2))
    peaks, props = find_peaks(avg_subcube_spectrum, height=avg_subcube_spectrum.max().value/4, width=0.1)

    if len(peaks)>1: #for multi-peaks, just fit the top two
        print("Double Peak")
        A0 = avg_subcube_spectrum[peaks][np.argsort(-props['peak_heights'])[0:2]]
        #mu0 = cube.spectral_axis[peaks][np.argsort(-props['peak_heights'])[0:2]] 
        #sig0 = props['widths'][np.argsort(-props['peak_heights'])[0:2]]
        
        ###instead, make both of the m0 = centroid of EVF from table, same with sig0 using vel extent from table!
        mu0 = subcube.spectral_axis[peaks][np.argsort(-props['peak_heights'])[0:2]] 
        sig0 = props['widths'][np.argsort(-props['peak_heights'])[0:2]]
        

        params0 = [mu0,A0,sig0]
        print(params0)
        params,cov=curve_fit(double_gauss,subcube.spectral_axis,avg_subcube_spectrum,p0=[mu0.value[0],mu0.value[1],A0[0].value,A0[1].value,sig0[0],sig0[1]])

        mu1,mu2,A1,A2,sig1,sig2 = params
        print('mu 1 = ',mu1,'mu 2 = ', mu2)
        print('sigma 1 = ',sig1,'sigma 2 = ', sig2)

        fig, ax = plt.subplots(1, 1)
        ax.plot(subcube.spectral_axis,avg_subcube_spectrum)
        #ax.plot(subcube.spectral_axis[peaks][np.argsort(-props['peak_heights'])[0:2]], avg_subcube_spectrum[peaks][np.argsort(-props['peak_heights'])[0:2]], "x")
        ax.plot(subcube.spectral_axis,double_gauss(subcube.spectral_axis.value,mu1,mu2,A1,A2,sig1,sig2),'--', c='r')
        ax.text(.25, .95, 'EVF ID: {}'.format(ID_num), fontsize=13, ha='right', va='top', 
            transform=ax.transAxes)
        
        ax.text(.95, .95, 'FWHM 1: {:.2f} km/s'.format(abs(sig1*2.355)), fontsize=13, ha='right', va='top', 
            transform=ax.transAxes)
        ax.text(.95, .85, 'FWHM 2: {:.2f} km/s'.format(abs(sig2*2.355)), fontsize=13, ha='right', va='top', 
            transform=ax.transAxes)
        ax.set_xlabel('km/s')
        ax.set_ylabel('Jy/beam')
        plt.savefig('./spectra/EVF_{}_{}_1Dspectrum_linewidths.png'.format(ID_num, EVF_filename), bbox_inches='tight')
        plt.show()
        return mu1,mu2,A1,A2,sig1,sig2
 
    
    else: #for single peak, one gaussian is enough
        print("Single Peak")
        A0 = max(avg_subcube_spectrum)
        #mu0 = cube.spectral_axis[peaks][np.argsort(-props['peak_heights'])[0]]
        #sig0 = props['widths'][np.argsort(-props['peak_heights'])[0]]
        ###instead, make both of the m0 = centroid of EVF from table, same with sig0 using vel extent from table!
        mu0 = EVF_v_cent[tab_ind]
        sig0 = EVF_v_width[tab_ind]
        
        params0 = [mu0,A0,sig0]
        params,cov=curve_fit(single_gauss,subcube.spectral_axis,avg_subcube_spectrum,p0=[mu0,sig0, A0.value])
        mu1,sig1,A1 = params

        print('mu 0 = ',mu0, 'mu 1 = ',mu1)
        print('A 0 = ',A0, 'A 1 = ',A1)
        print('sigma 0 = ',sig0, 'sigma 1 = ',sig1)
        fig, ax = plt.subplots(1, 1)
        ax.plot(subcube.spectral_axis,avg_subcube_spectrum)
        #ax.plot(cube.spectral_axis[peaks][np.argsort(-props['peak_heights'])[0]], avg_spectrum[peaks][np.argsort(-props['peak_heights'])[0]], "x")
        ax.plot(subcube.spectral_axis,single_gauss(subcube.spectral_axis.value,mu1,sig1,A1),'--', c='r')
        ax.text(.25, .95, 'EVF ID: {}'.format(ID_num), fontsize=13, ha='right', va='top', 
            transform=ax.transAxes)
        
        ax.text(.95, .95, 'FWHM 1: {:.2f} km/s'.format(abs(sig1*2.355)), fontsize=13, ha='right', va='top', 
            transform=ax.transAxes)
        ax.set_xlabel('km/s')
        ax.set_ylabel('Jy/beam')
        plt.savefig('./spectra/EVF_{}_{}_1Dspectrum_linewidths.png'.format(ID_num, EVF_filename), bbox_inches='tight')

        plt.show()
        return mu1, np.nan, A1, np.nan, sig1, np.nan  




EVF_ID_list, sigma1_list, sigma2_list, FWHM1_list, FWHM2_list = [], [], [], [], []

for file in glob.glob(CS_cubes_path + '*.fits', recursive=True):
    CS_cube = file
    ID_num = os.path.basename(file).split('_',2)[1]
    EVF_filename = os.path.basename(file).split('_',2)[2].split('.fits',1)[0]

    #read in the CS cube
    data = pyfits.open(CS_cube) 
    cube = SpectralCube.read(data)      
    data.close()      
    cube.allow_huge_operations=True
    mu1,mu2,A1,A2,sig1,sig2 = fit_gaussians(cube, ID_num,EVF_filename)

    if np.isnan(sig2) == False:
        FWHM_1, FWHM_2 = 2.355 * sig1, 2.355 * sig2
        sigma1_list.append(abs(np.round(sig1,3)))
        sigma2_list.append(abs(np.round(sig2,3)))
        FWHM1_list.append(abs(np.round(FWHM_1,3)))
        FWHM2_list.append(abs(np.round(FWHM_2,3)))
    else:
        FWHM_1 = 2.355 * sig1
        sigma1_list.append(abs(np.round(sig1,3)))
        sigma2_list.append(str('-'))
        FWHM1_list.append(abs(np.round(FWHM_1,3)))
        FWHM2_list.append(str('-'))
        
    #Make the CS PV diagrams
    make_pv_plots(CS_cube, ID_num,EVF_filename)
    

#Write output to table
EVF_linewidth_tab = Table()
EVF_linewidth_tab['EVF_ID'] = EVF_ID_list
EVF_linewidth_tab['sigma_1']   = sigma1_list
EVF_linewidth_tab['sigma_2']   = sigma2_list
EVF_linewidth_tab['FWHM_1']  = FWHM1_list
EVF_linewidth_tab['FWHM_2']  = FWHM2_list

EVF_linewidth_tab.write(basepath + 'spectra/EVF_linewidth_tab.tex', format = 'latex', overwrite=True)


##Now make PVs for the other lines##
for file in glob.glob(SiO_cubes_path + '*.fits', recursive=True):
    SiO_cube = file
    ID_num = os.path.basename(file).split('_',2)[1]
    EVF_filename = os.path.basename(file).split('_',2)[2].split('.fits',1)[0]

    #read in the SiO cube
    data = pyfits.open(SiO_cube) 
    cube = SpectralCube.read(data)      
    data.close()      
    cube.allow_huge_operations=True   
    make_pv_plots(cube, ID_num,EVF_filename)

for file in glob.glob(HC3N_cubes_path + '*.fits', recursive=True):
    HC3N_cube = file
    ID_num = os.path.basename(file).split('_',2)[1]
    EVF_filename = os.path.basename(file).split('_',2)[2].split('.fits',1)[0]

    #read in the HC3N cube
    data = pyfits.open(HC3N_cube) 
    cube = SpectralCube.read(data)      
    data.close()      
    cube.allow_huge_operations=True   
    make_pv_plots(cube, ID_num,EVF_filename)