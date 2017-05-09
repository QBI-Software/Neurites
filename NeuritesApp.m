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
% Version: 2.0 (Oct 2015) - Synapse analysis
% Version: 3.0 (Aug 2016) - Annulus analysis
% Developer: Liz Cooper-Williams,e.cooperwilliams@uq.edu.au
% Copyright {2017} { Queensland Brain Institute, University of Qld }
% 
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
% 
%        http://www.apache.org/licenses/LICENSE-2.0
% 
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function NeuritesApp()

  clear
  close all
  format compact
  feature('jit',0)
  feature('accel',0)
  %add paths
  addpath(pwd)
  subd = fullfile(pwd, 'methods');
  addpath(subd)
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