%%Test NeuritesAnalyser and NeuritesSynapse
imgname='100715 annotated.tif';%example.tif'
cell1file = 
I = imread(imgname);
% thresholds – 1x3 matrix with the threshold values
% for the R, G and B color channels - from histogram
[cell1,cell2] = SplitCells(I,256);
N = NeuritesAnalyser(imgname,'tmp_roi.tif',cell1,cell2);
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

N = N.measureSynapses(showtypes,cell1file,cell2file);
