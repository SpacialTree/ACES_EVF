# ACES_EVF_HVCC
This repository holds the analysis code used to search for and analyze High Velocity Compact Cloud candidates in the ACES CS (2-1) data.

## products 
The aces_evf/products/ directory contains the following files:
**'HVCC_regions.reg' :** ds9 format region file with spatial extents of HVCCs from the catalog
**'physical_properties_paper_tab_v5.tex' :** Table 2 from the paper (has not been fully converted to machine readable format yet)
**High resolution cutout cubes of the HVCCs are available on Globus**

## Other files
The code used for analysis is in the main aces_evf/ directory, including the following files:

**'ACES_CS_HVCC_v2.ipynb'** : Dendrogram extraction of HVCCs in the ACES CS (2-1) cube	
**'ACES_CS_HVCC_lbplot.ipynb'** : Figure 2	
**'CS_HVCC_catalog.py'** : 	Generation of catalog from the CS dendrogram  **'CS_dendro.py'** was the original test version.
**'CS_dendro_match.py'** : cross matching dendrogram detections and by-eye detections		

**'line_ratio_histograms.ipynb'** : stack histograms of line ratios colored by PV morphological class
				
**'make_mom0_and_ratio-maps_3sigma_masking.py'**, **'make_mom0_and_ratio-maps_3sigma_masking_HNCO_MOM0_ONLY.py'** : Creation of mom0 and line ratio maps for all ACES lines **'make_mom0_and_ratio_maps.py'** was the original test version before breaking up into the 3sigma masks.

**'50kmsC_HVCCs.ipynb'** : Exploring overlap of HVCCs on the 50 km/s cloud and CND area in spatial and velocity space. Used for creating Figures 9 and 10 from the paper. 				

**'make_pv_linewidth.ipynb'** : calculation of FWHM values and statistics (Figs 5, 7, 11, 12) 

**'pv_maker.py'** : Creation of PV slices for visual inspection

**'PPV_visual.ipynb'** : Initial PPV visualization. See updated webapp here: https://github.com/danilipman?tab=repositories 						

**'oka_comp_table.ipynb'** : Comparison with Oka 2012 and Oka 2022 tables

**'evfs_selected.ipynb'** : Used to create: average CS(2-1) spectrum, Position-Velocity Diagram in Galactic latitude and Galactic longitude, and Moment 0 map for selected EVFs.

## Other misc code used during paper sprint : 
Histogram of Line Ratios PV 
Classifications-Final Version.ipynb	products
pv_maker.ipynb
cutout_pv.ipynb							
evf_analysis.ipynb						
pv_plotter.ipynb
evf_cutouts.ipynb						
pv_plotter.py
evf_cutouts.py							
region_list.reg
explore_pv.ipynb						
rename_resort.py
inflow_EVF_PPV_maker.ipynb
make_pv_linewidth.py
oka_table.ipynb

