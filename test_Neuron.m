%Test for loadCSVtree
csvfile = 'D:\Projects\Williams_neurites\Alpha Ganglion cells (Arne)\annulusdata\16_03_08cell01_a.csv';
%rgbimg = 'D:\Projects\Williams_neurites\Alpha Ganglion cells (Arne)\annulusdata\16_03_08 cell 01 0.1 mum per pixel.tif';
%I = imread(rgbimg);
%n = NeuritesAnnulusAnalyser(I);
%soma = n.soma
configfile = 'D:\Projects\Williams_neurites\Alpha Ganglion cells (Arne)\annulusdata\neurites_annulus_config.csv';
config = readtable(configfile)
soma = NeuritesSoma(0,0,[config.CentroidX config.CentroidY])
neuron = loadCSVTree(csvfile,soma,config.Scale);