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
olist = {};
sprintf('%d levels',length(branches))
out = [];
out2 = [];
out3 = [];
N = branches{1,1}(1,1); %rootnode
out2 = inOrder(data, N, out2);
out = postOrder(data, N, out); %THIS ONE BEST
out3 = preOrder(data, N, out3);

%branchdistances = arrayfun(@(n) n.branchlength, out)
%bcount = length(branchdistances);
%bpair = [1 2;3 4;11 5;12 6;7 8;13 14;10 15;16 9] - this works but HOW??

orderedout = arrayfun(@(n) n.isBranchpoint(),out)
bps = out(orderedout)
eps = out(~orderedout)

branchdistances = cat(1, (arrayfun(@(n) n.branchlength,eps)),(arrayfun(@(n) n.branchlength,bps)))
branchleaves = arrayfun(@(n) n.nodelevel,eps)
branchnames = arrayfun(@(n) sprintf('Branch %d',n.nodelevel),bps, 'uniformoutput', false)
branchdata = getBranch2Leaf(out)
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
      for m=1:length(branches{1,j})
        N = branches{1,j}(m,1);
        data{N.id} = N;
      end
    end
end

function output = inOrder(data, nnode, output)
    if (~isempty(nnode.leftnode))
        output = inOrder(data, data{nnode.leftnode.id}, output);
    end
    output = cat(1, output, nnode);
    if (~isempty(nnode.rightnode))
        output = inOrder(data, data{nnode.rightnode.id}, output);
    end
end

function output = postOrder(data, nnode, output)
    if (~isempty(nnode.leftnode))
        output = postOrder(data, data{nnode.leftnode.id}, output);
    end
    
    if (~isempty(nnode.rightnode))
        output = postOrder(data, data{nnode.rightnode.id}, output);
    end
    output = cat(1, output, nnode);
end
    
function output = preOrder(data, nnode, output)
    output = cat(1, output, nnode);
    if (~isempty(nnode.leftnode))
        output = postOrder(data, data{nnode.leftnode.id}, output);
    end
    
    if (~isempty(nnode.rightnode))
        output = postOrder(data, data{nnode.rightnode.id}, output);
    end
    
    
end

function bdata = getBranch2Leaf(nlist)
    bdata = [];
    for i=1:length(nlist)
        leaf =i;
        if (nlist(i).hasParent())
            pid = nlist(i).parentnode.id;
            parent = arrayfun(@(x) ~isempty(find(x.id == pid)), nlist,'uniformoutput',true)
            branch = find(parent) %should give idx
            bdata = cat(1,bdata,[branch leaf])
        end
    end
end
