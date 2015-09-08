QBI Williams Neurites Synapse Analyser
---------------------------------------
Description: Developed for Simon de Croft (Williams lab) for the 
detection and analysis of synaptic regions of two neurons.

Input: 
-----
   1. Image (tiff or jpg) with both neurons as two color tracing (from Neurolucida)
   2. CSV analysis files for each cell 
      ie., Tab 2 'Segment Points to Dendrites') labelled accordingly 
      eg. DSdata.csv and SBACdata.csv

Output:
-----
  1. Neurites_data.csv (analysis)
  2. Neurites_roi.tif (ROI mask image)
  3. Neurites_config.csv (Configuration for overlay of CSV data)

Requires:
-----
  Matlab Image Processing Toolbox
  Matlab Statistics Toolbox

Usage:
-----
  1. Load image file (tif or jpg)
  2. Load analysis files (csv) as above
  3. Set ROI with ROI tool (radio button in Actions->Image view) -
      double-click to save
  4. Register CSV overlay - adjust with Configuration options then 
  click on 'Save config' this will reload when analysis files are loaded:
    Scale: integer (eg 1, 10)
    Fit tolerance: Number of pixels either side - usually 1
    ShiftX: Number of pixels to shift along x-axis (can be
          determined with 'Image Tool -> pixel region button or
          measurement tool')
    ShiftY: Number of pixels to shift along y-axis (as for x)
    Min length: Minimum length of a synaptic region (eg 10 pixels)
    Tolerance: Nearest neighbour matching of apposing regions (eg 1
          pixel)
  5. Identify Synaptic Regions - runs analysis. Check results in windows and table.

NB. Output files will be overwritten on each run - to save just copy file
directly

