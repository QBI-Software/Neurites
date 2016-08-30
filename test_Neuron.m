function test_Neuron()
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
sprintf('Loaded neuron with %d trees', length(neuron));
%INTEGRATE: lookup from analyse app via id == fR
%TODO: Create phytree
branches = neuron{1,1};
bcount = length(branches);
sprintf('Num branches: %d',length(branches)) 
sprintf('Num nodes: %d',bcount)
%load tree vars
% branchdata: contains branch pairs by IDs eg[1 2;3 4]; - levels
% branchdistances: 
% M =
% 
%   (10,1)        1
%   (10,2)        1
%   (11,3)        1
%   (11,4)        1
%   (12,5)        1
%   (13,6)        1
%   (14,7)        1
%   (14,8)        1
%   (17,9)        1
%   (16,10)       1
%   (12,11)       1
%   (13,12)       1
%   (15,13)       1
%   (15,14)       1
%   (16,15)       1
%   (17,16)       1
% 
% 
% ID = 
% 
%     '4A'
%     '4B'
%     '7A'
%     '7B'
%     '6'
%     '5'
%     '5A'
%     '5B'
%     '2'
%     'Branch 1'
%     'Branch 2'
%     'Branch 3'
%     'Branch 4'
%     'Branch 5'
%     'Branch 6'
%     'Branch 7'
%     'Branch 8'
% 
% 
% DIST =
% 
%   123.3000
%   133.1000
%   125.0000
%   290.7000
%   324.2000
%   371.3000
%   323.0000
%   184.5000
%   178.5000
%   159.4000
%   139.7000
%   280.9000
%    23.0000
%    89.9000
%    31.6000
%    24.9000
%    50.0000
% br = get(fulltree) -> load n0001.tree and export to workplace as fulltree
%Loop through each level - find branchpoints then endpoints
data = loadById(branches);
bcount = 0;
branchptrs = [];
branchleaves = [];
branchnames = [];
branchdistances = [];% [290.7;125;324.2;139.7;0]
%generate list of endpoints
eps = [];
bps = [];
ptrs = [];
%for j=[length(branches):-1:1]
for j=1:length(branches)
    %Find end branches - first string
    for (m=1:length(branches{1,j}))
        N = branches{1,j}(m,1)
        if (N.isBranchpoint())
            bps = cat(1,bps,N);
            ptrs = cat(1, ptrs,[N.countid N.leftnode.countid]);
            ptrs = cat(1, ptrs,[N.countid N.rightnode.countid]);
        else
            eps = cat(1,eps,N);
        end
    end
  
end
branchdistances = cat(1, (arrayfun(@(n) n.branchlength,bps)),(arrayfun(@(n) n.branchlength,eps)))
branchleaves = arrayfun(@(n) n.nodelevel,eps)
branchnames = arrayfun(@(n) sprintf('Branch %d',n.nodelevel),bps, 'uniformoutput', false)
tree = phytree(branchdata,branchdistances)
names = get(tree,'LeafNames')

%Show plot
view(tree);
newick = getnewickstr(tree); %save as newick format for other programs
%h = plot(tree,'orient','top');
xlabel('Distance from soma (um)')
%set(h.terminalNodeLabels,'Rotation',65)
end

function data = loadById(branches)
    data = {};
    for j=1:length(branches)
      for (m=1:length(branches{1,j}):1)
        N = branches{1,j}(m,1);
        data{N.id} = N;
      end
    end
end