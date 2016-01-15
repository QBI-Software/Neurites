%%Test NeuritesAnalyser and NeuritesSynapse
imgname='Neurons.tif';%'100715 annotated.tif';%example.tif'
imgroi ='neurites_roi.tif';
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
cell1file='DSdata.csv';
cell2file='SBACdata.csv';
scale=2.84;
shiftx=0;
shifty=-935;
fit=10;
cell1label='DS';
cell2label='SBAC';
N = N.measureSynapses(showtypes,cell1file,cell2file,scale,shiftx,shifty,fit);
[colnames,T1] = N.generateCentredDataTable(showtypes,cell1label, cell2label);
%save to file
outputfile = fullfile(pathname, 'neurites_data_centred.csv');
saveDataFile(outputfile, colnames, T1);
