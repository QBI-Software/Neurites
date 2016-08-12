% QBI Williams Neurites Synapse Analyser
% 
% Description: Developed for the detection and analysis of dendritic fields.
%   Select: 1. Synapse Analysis :
%              Detects and analyses regions of apposition (potential synaptic
%              regions) of overlapping dendritic fields of two neurons.
%           2. Annulus Analysis :
%              Analyses the neurites of a single neuron under an annulus. 
%              Also presents information on a dendrogram with higlighted branches.
% 
% Requires:
%   Matlab Image Processing Toolbox
%   Matlab Statistics Toolbox
% 
% Code and further information is available at : 
%       https://github.com/QBI-Software/Neurites
% Please see Instructions for Usage at :
%       https://github.com/QBI-Software/Neurites/wiki/Instructions
% 
% Developed by Liz Cooper-Williams, QBI (e.cooperwilliams@uq.edu.au)
% Version: 2.0 (Oct 2015) - Synapse analysis
% Version: 3.0 (Aug 2016) - Annulus analysis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NeuritesApp()

  clear
  close all
  format compact
  feature('jit',0)
  feature('accel',0)
  %add paths
  addpath(pwd)
  addpath('./methods')
  % Select application to run
  choice = questdlg('Select Application to run', ...
	'Neurites analysis app', ...
	'Synapse analysis','Annulus analysis','Cancel','Synapse analysis');
    % Handle response
    switch choice
        case 'Synapse analysis'
            disp([choice ' loading.'])
            %Launch UI
            NeuritesAppUI()
        case 'Annulus analysis'
            disp([choice ' loading.'])
            %Launch UI
            NeuritesAnnulusAppUI()
        case 'Cancel'
            disp('No app selected.')
    end
  
  
 

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end