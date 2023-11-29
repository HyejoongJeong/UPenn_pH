# UPenn_pH

MICROWELL ARRAY DROPLET IMAGE ANALYSIS

This repository contains ImageJ and MATLAB files used in the manuscript {INSERT FINAL MANUSCRIPT TITLE}. The files are listed in the order they must be run.

imageProcessingMacro contains code for automated particle tracking of cell bodies within microwell array images.

Droplet_analysis_driverCode contains driver code to automate image analysis across multiple images. This code calls fxn_droplet_analysis using multi-channel microwell array images over a time course as input and as output provides single, empty, and multi cell droplet indices as well as cell and droplet fluorescence intensities over time.

fxn_droplet_analysis contains code for droplet segmentation, cell counting and tracking, and fluorescence signal extration.

Droplet_analysis_generatePlots contains code to generate various plots of fluorescence intensity data. Plots include GFP intensity of cells and droplets over time and Cy5/Cy3 intensity ratios (used in pH calculation) in droplets over time. These data are grouped by droplet classifications of empty, single-, and multi- cell droplets/wells. 

Evaluate_all_positions contains code to combine and generate plots for multiple images within the same data set. 
