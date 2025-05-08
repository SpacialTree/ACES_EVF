#!/usr/bin/env python
# coding: utf-8

# In[3]:


#from IPython.display import display, HTML
#display(HTML("<style>.container { width:96% !important; }</style>"))


# In[2]:


import os, glob
from pathlib import Path

import numpy as np
import pandas as pd
import pylab
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table


import matplotlib.pyplot as pl
import matplotlib.colors as mc
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from mpl_toolkits.mplot3d import axes3d, art3d, proj3d
from matplotlib import colors
from matplotlib.colors import Normalize
import matplotlib.cm as cm


from scipy.interpolate import NearestNDInterpolator
from numpy import linspace, array, logspace, sin, cos, pi, arange, sqrt, arctan2, arccos
from adjustText import adjust_text
import matplotlib.patheffects as PathEffects
from astropy.coordinates import Angle, SkyCoord, Longitude
from astropy.nddata.utils import Cutout2D
import astropy.units as u
from astropy.wcs.utils import pixel_to_skycoord
import astropy.io.fits as pyfits


from regions import Regions
from spectral_cube import SpectralCube
from spectral_cube import Projection

from reproject import reproject_interp




#Paths and definitions
drivepath = '/orange/adamginsburg/ACES/mosaics/cubes/' #Data directory containing cubes
mom0dir = '/orange/adamginsburg/ACES/broadline_sources/EVFs/moments/'#Directory where mom0 maps are stored. Please don't put other things in this directory.
savedir_figure = '/orange/adamginsburg/ACES/broadline_sources/EVFs/ratios/' #Directory where ratio map figures will be saved (with folders in this dir made for each EVF)


EVF_tab = Table.read('/blue/adamginsburg/savannahgramze/ACES_EVF/aces_evf/Filtered_EVFs_table.ecsv')
EVF_reg = Regions.read('/blue/adamginsburg/savannahgramze/ACES_EVF/aces_evf/EVF_reg_list.reg')




vel_range_list = []
delta_l_list  = []
delta_b_list = []
lb_list = []
for evf in EVF_tab:
    vel_range = (evf['min_v'], evf['max_v'])
    vel_range_list.append(vel_range)
    delta_l_list.append(evf['deltal'])
    delta_b_list.append(evf['deltab'])
    lb=(evf['l'], evf['b'])
    lb_list.append(lb)

linetracers = ['CS21', 'H13CN', 'H13COp', 'SiO21', 'SO32', 'SO21', 'HN13C', 'HC3N', 'HNCO_7m12mTP']#, 'HCOP_noTP']


#Make and save moment 0 maps for all line tracers
for line in linetracers:
    mom0folder = mom0dir + '{}'.format(line)
    Path(mom0folder).mkdir(parents=True, exist_ok=True) #look for folder for the linetracer mom0 maps, if it doesn't exit, make it
    for file in glob.glob(drivepath + line + '_CubeMosaic.fits', recursive=True):
        data = pyfits.open(file) 
        cube = SpectralCube.read(data)
        data.close()
        cube.allow_huge_operations=True

        for i in range(len(EVF_reg)):
            subcube = cube.subcube_from_regions([EVF_reg[i]])
            subcube = subcube.spectral_slab(vel_range_list[i][0]* u.km / u.s, vel_range_list[i][1]* u.km / u.s)
            mom0 = subcube.moment(order=0)


            mom0.write(mom0folder +f'/l{lb_list[i][0]}_b{lb_list[i][1]}_{line}mom0.fits', overwrite=True)


#Gets an array of all the EVF sources to make ratio maps for
file_names = []

for line in linetracers:
    mom0folder = mom0dir + '{}'.format(line)
    for file in os.listdir(mom0folder):
        filename = os.fsdecode(file)
        if filename.endswith('.fits'):
            file_names.append(filename)
            
#Gets an array of just the EVF name strings (for labels/file naming/matching)
source_names=[]
i=0
while i<len(file_names):
    array = file_names[i].split('_')
    if len(array)==3:
        l = array[0] #This is how I piece together the EVF source name from file name
        b = array[1]
        source = l + '_' + b
        source_names.append(source)
    else: 
        print("There's a file name that's structured differently than others in this directory!")
    i+=1

sources_listed_once = [i for n, i in enumerate(source_names) if i not in source_names[:n]] #removes any duplicates from list (& keeps order)
print("The number of EVFs to do ratio maps for: ", len(sources_listed_once))



#This will make and save line tracer ratio maps for all EVF sources
noise_threshold = 9*10**(-3) #Jy/beam noise threshold below which we mask in our ratio maps (including any absorption)
s=0
i=0
while s<len(sources_listed_once): #Iterates through each source
    #Figure out which lines are available for each EVF source
    source=sources_listed_once[s]
    lines = linetracers.copy() #array for line names
    mom0paths = [] #array for the source's mom0 maps
    for line in lines: 
        for folder in os.listdir(mom0dir): #goes through line folders (and any other files/folders) in mom0dir
            foldername = os.fsdecode(folder)
            if line in foldername:
                print(line)
                for file in os.listdir(mom0dir+foldername):
                    filename = os.fsdecode(file)
                    if source in filename: 
                        print(filename)
                        mom0paths.append(filename)
    
    lines_avoidredundant = lines.copy() #arrays I'll remove elements from to avoid redundant pairs
    mom0_avoidredundant = mom0paths.copy()
    #print("Source: ")
    #print(sources_listed_once[s])
    
    i=0
    while i<len(lines):
        r=0
        while r<len(lines_avoidredundant):
            if lines[i] == lines_avoidredundant[r]: #don't bother calculating CS-CS, HCN-HCN, etc, ratios
                r+=1
            else: #run everything else
                print(lines[i], "-->", lines_avoidredundant[r])
                #Open mom0 map of line 1:
                mom0path1 = mom0dir + lines[i] + '/' + mom0paths[i]
                #print(mom0path1)
                hdul = fits.open(mom0path1)
                sc_line1_moment0 = Projection.from_hdu(hdul) #makes 2D spectral cube object called a "Projection"
                
                #Open mom0 map of line2:
                mom0path2 = mom0dir + lines_avoidredundant[r] + '/' + mom0_avoidredundant[r]
                #print(mom0path2)
                hdul = fits.open(mom0path2)
                sc_line2_moment0 = Projection.from_hdu(hdul) #makes 2D spectral cube object called a "Projection"
                
                #Match shapes of 2D moment maps so we can calculate ratios
                sc_line1_moment0_reproject, footprint = reproject_interp(sc_line1_moment0.hdu,sc_line2_moment0.header) #reproject line1 to line2

                #Now that images are same size, compute ratio of moment maps
                ratio = sc_line2_moment0.hdu.data/sc_line1_moment0_reproject

                #Mask pixels (exclude certain values) in the ratio map:
                badpix = pylab.where(sc_line1_moment0_reproject<noise_threshold) # Identify emission below a threshold to mask for line 1
                badpix2 = pylab.where(sc_line2_moment0.hdu.data<noise_threshold)   #Line2 noise and absorption
                ratio[badpix] = np.nan # Mask the ratio map
                ratio[badpix2] = np.nan

                
                #Plot (with a color bar)
                if np.isnan(np.nanmin(ratio))==True or np.isnan(np.nanmax(ratio))==True:
                    print("!!!!! ALL-NAN RATIO MAP !!!!!!!!!!!!")
                    print("Check mom0 maps for: ", source)
                    print("Potentially problematic maps: ", mom0path1, mom0path2)
                    
                else: 
                    normalizer=colors.LogNorm(vmin=np.nanmin(ratio),vmax=np.nanmax(ratio)) #log scaling on ratio map and color bar
                    imnorm=cm.ScalarMappable(norm=normalizer, cmap='rainbow')

                    fig1 = pylab.figure(1,figsize=(15,15)) #Figure size

                    ax1 = pylab.subplot(1,1,1,projection=sc_line2_moment0.wcs) 
                    im1 = pylab.imshow(ratio,cmap='rainbow', norm=normalizer)


                    gallong = ax1.coords[0]                                                                
                    gallat = ax1.coords[1]
                    gallong.set_ticks(size=-3) #Axis ticks                                                                                      
                    gallat.set_ticks(size=-3)  

                    pylab.title('ACES %s' %source + ' Line Ratio Map', fontsize=20) 
                    pylab.xlabel('Galactic Longitide',fontsize=15,labelpad=1)                               
                    pylab.ylabel('Galactic Latitude',fontsize=15,labelpad=0)
                    ax1.tick_params(axis = 'both', which = 'major', labelsize = 15)
                    #Colorbar
                    cb=pylab.colorbar(imnorm,fraction=0.046,pad=0.04)                                      
                    cb.set_label(label='%s' %lines_avoidredundant[r] + ' /%s Ratio' %lines[i], fontsize=15,rotation=270,labelpad=20)
                    cb.ax.tick_params(which = 'major', labelsize = 10)   
                    #Save each figure
                    Path(savedir_figure+source).mkdir(parents=True, exist_ok=True) ##look for folder individual EVF, if it doesn't exit, make it
                    savepath_figure = savedir_figure + source + '/RatioMap_' + lines_avoidredundant[r] + '_' + lines[i] + '_' + source + '.png'
                    pylab.savefig(savepath_figure)
                r+=1
        del lines_avoidredundant[0] #remove 1st element each line iteration to redundant runs (if CS-HC3N done, don't do HC3N-CS)
        del mom0_avoidredundant[0]
        i+=1
    #print("**************************************")

    s+=1
