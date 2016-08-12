QBI Williams Neurites Synapse Analyser
---------------------------------------
Description: Developed for the Williams lab for the morphological
detection and analysis of dendritic fields.
	1. Synapse Analysis :
       Detects and analyses regions of apposition (potential synaptic        regions) of overlapping dendritic fields of two neurons.
	2. Annulus Analysis :
       Analyses the neurites of a single neuron under an annulus. 
       Also presents information on a dendrogram with higlighted branches.

Input: 
-----
   1. Image (tiff or jpg) (eg, from Neurolucida tracing)
   		a. Synapse Analysis: requires two neurons in separate colors 
   		b. Annulus Analysis: single neuron
   2. CSV analysis files for each cell 
      eg, Neurolucida measurements: Tab 2 'Segment Points to Dendrites' 

Output files:
NB. Output files will be overwritten on each run (rename files in browser to save if necessary)

Synapse Analysis Output:
-----
  1. Neurites_data.csv (analysis)
  2. Neurites_roi.tif (ROI mask image)
  3. Neurites_config.csv (Configuration for overlay of CSV data)
  
Annulus Analysis Output:
-----
  1. Neurites_data.csv (analysis)
  2. Neurites_roi.tif (ROI mask image)
  3. Neurites_config.csv (Configuration for overlay of CSV data)

Requires:
-----
  Matlab version R2015a (other versions not tested)
  Matlab Image Processing Toolbox
  Matlab Statistics Toolbox

Usage:
-----
  Refer to https://github.com/QBI-Software/Neurites/wiki/Instructions



