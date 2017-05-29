%%Test NeuritesAnalyser and NeuritesSynapse
%GENERATE IMAGE WITH CSV data ONLY
cell1file='sampledata/010415_DSdata.csv';
cell2file='sampledata/010415_SBACdata.csv';
samplecells = {cell1file, cell2file};
for m=1:length(samplecells)
    cellfile = samplecells{1,m};
    T = readtable(cellfile);
    centroidxy = [mean(T.StartX) mean(T.StartY)]; %estimate from CSV data
    scale = 1;
    soma = NeuritesSoma(0,0,centroidxy);
    tree = loadCSVTree(cellfile, soma, scale);
    sprintf('Loaded neuron from %s with %d trees', cellfile,length(tree));
    %Actual centroid data - reload tree
    [cx,cy] = findcentroid(tree);
    soma = NeuritesSoma(0,0,[cx,cy]);
    tree = loadCSVTree(cellfile, soma, scale);
    samplecells{2,m} = tree;
    plotdigineuron(tree);
end
    
% T1 = readtable(cell1file);
% T2 = readtable(cell2file);
% %config = readtable('sampledata/neurites_config.csv');
% centroidxy = [mean(T1.StartX) mean(T1.StartY)]; %estimate from CSV data
% scale = 1;
% soma1 = NeuritesSoma(0,0,centroidxy);
% tree1 = loadCSVTree(cell1file, soma1, scale);
% sprintf('Loaded neuron 1 with %d trees', length(tree1));
% centroidxy = [mean(T2.StartX) mean(T2.StartY)]; %estimate from CSV data
% soma2 = NeuritesSoma(0,0,centroidxy);
% tree2 = loadCSVTree(cell2file, soma2, scale);
% sprintf('Loaded neuron 2 with %d trees', length(tree2));
% 
% %Revise centroid data
% [cx,cy] = findcentroid(tree1);
% soma1 = NeuritesSoma(0,0,[cx,cy]);
% tree1 = loadCSVTree(cell1file, soma1, scale);
% [cx,cy] = findcentroid(tree2);
% soma2 = NeuritesSoma(0,0,[cx,cy]);
% tree2 = loadCSVTree(cell2file, soma2, scale);
% plotdigineuron(tree1);
% plotdigineuron(tree2);
%%TODO: DIGITAL NEURON: loop through other tree until find intersect C = intersect(A,B,'rows')


%%%%METHODS FOR NEURON TREE%%%%%%
function [cx,cy] = findcentroid(tree1)
    polycoordsX = [];
    polycoordsY = [];
    for k=1:length(tree1)
        N = tree1{k,1}{1,1}; %root node
        if (N.nodelevel == 1)
            polycoordsX(end+1) =N.points(1,1);
            polycoordsY(end+1) =N.points(1,2);
            
        end
    end
    %if lines don't cross - use polycoords
    %[xi,yi] = polyxpoly(x1,y1,x2,y2) returns the intersection points of two polylines in a planar, Cartesian system
    n = ceil(length(polycoordsX)/2);
    x = polycoordsX;
    y = polycoordsY;
    %Append first set of point to complete polygon
    x(end+1) = x(1);
    y(end+1) = y(1);
    x1 = x(1:n);
    y1 = y(1:n);
    x2 = x(n+1:end);
    y2 = y(n+1:end);
    [xi,yi] = polyxpoly(x1,y1,x2,y2);
    if(length(xi) > 1) %should at least match starting point
    	
        %convert to polar to sort coords
        [th,r]=cart2pol(polycoordsX,polycoordsY);
        %sort by angle
        mapObj = containers.Map(th,r)
        %polarplot(cell2mat(mapObj.keys()),cell2mat(mapObj.values()))
        [x,y] = pol2cart(cell2mat(mapObj.keys()),cell2mat(mapObj.values()))
        %Append first set of point to complete polygon
        x(end+1) = x(1);
        y(end+1) = y(1);
    end
    
    plot(x,y,'LineStyle',':'); 
    cx = mean(x);
    cy = mean(y);
    hold on;
    plot(cx,cy,'marker','o','color','b');
end


function plotdigineuron(tree1)
    for k=1:length(tree1)

        branches1 = tree1{k,1};
        bcount = length(branches1);
        sprintf('Num branches: %d',length(branches1)) 
        sprintf('Num nodes: %d',bcount)
        hold on;
        data = loadById(branches1);
        N = branches1{1,1}(1,1); %root node
        out = [];
        out = postOrder(data, N, out); %THIS ONE BEST

        %%Plot branches - color coordinated levels
        colorstring = 'kbgrymckbgrymc';%colorstring(out(i).nodelevel)

        for i=1:length(out)
            plot(out(i).points(:,1),out(i).points(:,2),'color',colorstring(k)); 
        end
    end
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


    