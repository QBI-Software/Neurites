% QBI Williams Neurites Synapse Analyser
% 
% Aim: Developed for the detection and analysis of synaptic regions of 
%   two overlapping neurons.
% 
%   
% Input: 
%    1. Image (tiff or jpg) with both neurons as two color tracing (from Neurolucida)
%    2. CSV analysis files for each cell 
%       ie., Tab 2 'Segment Points to Dendrites') labelled accordingly 
%       eg. DSdata.csv and SBACdata.csv
% Description of analysis:
%    A single Neurolucida drawn image of two neurons is split based on color 
%   (detected by a single color per neuron only). The images are converted
%   to binary formats for analysis.  Boundary analysis from the Matlab 
%   Image Processing Toolbox is used to detect neurites within a manually
%   drawn ROI (saved as an image mask).  Each boundary segment is matched
%   against corresponding segments from the other cell and matching regions
%   (within a configurable tolerance) are detected as synaptic regions.
%   Data from Neurolucida measurements are used to populate metadata of the
%   synaptic regions.  The data (referred to as CSV data) is initially overlaid on
%   the image to obtain the scaling parameters.  Each synaptic
%   region is referenced to the CSV data to obtain distance and neurite
%   length measurements.  Duplicate synaptic regions are determined via the
%   range configured by the fit parameter.
%   Review of the detected synaptic regions displays the CSV data along
%   the image neurites.  The centroid of the somas of each cell is detected
%   via boundary and regionproperty analysis from the Matlab Image
%   Processing Toolbox and Matlab Statistics Toolbox.  The polar
%   coordinates of the medians of each cell's synaptic regions are
%   calculated relevant to that cell's soma centroid and mapped on a
%   compass map.
% Output:
%   1. Neurites_data.csv (analysis) and Neurites_data_review.csv (reviewed)
%   2. Neurites_roi.tif (ROI mask image)
%   3. Neurites_config.csv (Configuration for overlay of CSV data)
% 
% Requires:
%   Matlab Image Processing Toolbox
%   Matlab Statistics Toolbox
% 
% Usage: 1. Load image file (tif or jpg)
%       2. Load analysis files (csv) as above
%       3. Set ROI with ROI tool (radio button in Actions->Image view) -
%       double-click to save or load from previous run via File menu
%       4. Register CSV overlay: plots data from csv files onto image. This
%           needs to be fitted via following Configuration variables:
%           Scale: integer (eg 1, 10)
%           Fit tolerance: Number of pixels either side. 
%               Also used in detecting duplicate synapses (default 10)
%           ShiftX: Number of pixels to shift along x-axis (can be
%               determined with 'Image Tool -> pixel region button or
%               measurement tool')
%           ShiftY: Number of pixels to shift along y-axis (as for x)
%           Min length: Minimum length of a synaptic region (eg 10 pixels)
%           Tolerance: Nearest neighbour matching of apposing regions (eg 1
%           pixel)
%           - 'Save config' this will then reload automatically when CSV 
%           analysis files are loaded or load a file via File menu
%       5. Run analysis: Check results in windows and table.
%       6. Review: Allows examination of each synapse with paths to soma
%       and end of neurite.  Can remove erroneous branches and delete
%       duplicate synapses.  This is also good for printing.
%       7. Compass: Plots angles for each cell (direct line from synapse to
%       soma centroid)
%
% NB. Output files will be overwritten on each run - to save just copy file
% directly
% 
% Developed by Liz Cooper-Williams, QBI (e.cooperwilliams@uq.edu.au)
% Version: 2.0 (14 Oct 2015)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NeuritesApp()

  clear
  close all
  format compact
  feature('jit',0)
  feature('accel',0)
  %Launch UI
  NeuritesAppUI()
  
 

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end