%%Test NeuritesAnalyser and NeuritesSynapse
imgname='sampledata/010415 tracing A.tif';%'100715 annotated.tif';%example.tif'
imgroi ='sampledata/neurites_roi.tif';
I = imread(imgname);
% thresholds – 1x3 matrix with the threshold values
IG = rgb2gray(I);
% for the R, G and B color channels - from histogram
[cell1,cell2] = SplitCells(I,256);
centred = 1;
N = NeuritesAnalyser(imgname,imgroi,cell1,cell2);
minlength = 10;
tolerance = 2;
N = N.findSegments(minlength,tolerance);
figure
showplots = 1;
showtypes = [1 2 3];
N = N.analyseSynapses(showplots,showtypes);
ctr = 0;
for i=1:length(N.Synapses)
    syn = N.Synapses{i};
    if (syn.SynapseType == 1)
        ctr = ctr + 1;
    end
end
sprintf('Found %d synapses', ctr);
cell1file='sampledata/010415_DSdata.csv';
cell2file='sampledata/010415_SBACdata.csv';
%Load neurites_config.csv
M = readtable('sampledata/neurites_config.csv');

scale=M.Scale;
shiftx=M.Shiftx;
shifty=M.Shifty;
fit=M.Fit;
cell1label=M.Cell1;
cell2label=M.Cell2;

N = N.measureSynapses(showtypes,cell1file,cell2file,scale,shiftx,shifty,fit);
[colnames,T1] = N.generateCentredDataTable(showtypes,cell1label, cell2label);
%save to file
outputfile = fullfile('sampledata', 'neurites_data_centred.csv');
saveDataFile(outputfile, colnames, T1);
