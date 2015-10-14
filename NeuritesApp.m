% QBI Williams Neurites Synapse Analyser
% 
% Description: Developed for Simon de Croft (Williams lab) for the 
% detection and analysis of synaptic regions of two neurons.
% Input: 
%    1. Image (tiff or jpg) with both neurons as two color tracing (from Neurolucida)
%    2. CSV analysis files for each cell 
%       ie., Tab 2 'Segment Points to Dendrites') labelled accordingly 
%       eg. DSdata.csv and SBACdata.csv
%
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
  %Launch UI
  NeuritesAppUI()
  
 

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end