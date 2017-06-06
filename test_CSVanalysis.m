%%Test NeuritesAnalyser and NeuritesSynapse
%GENERATE IMAGE WITH CSV data ONLY
cell1file='sampledata/16-09-15_pair_04_AC.csv'; %010415_DSdata.csv';
cell2file='sampledata/16-09-15_pair_04_GC.csv'; %010415_SBACdata.csv';
%cell1file='sampledata/010415_DSdata.csv';
%cell2file='sampledata/010415_SBACdata.csv';
samplecells = {cell1file, cell2file};
colorstring = 'rgkbgrymckbgrymc';%colorstring(out(i).nodelevel)
for m=1:length(samplecells)
    cellfile = samplecells{1,m};
    T = readtable(cellfile);
    centroidxy = [mean(T.StartX) mean(T.StartY)]; %estimate from CSV data
    scale = 1;
    soma = NeuritesSoma(0,0,centroidxy);
    tree = loadCSVTree(cellfile, soma, scale);
    sprintf('Loaded neuron from %s with %d trees', cellfile,length(tree));
    %Actual centroid data - reload tree
    [cx,cy] = findsomacentroid(tree);
    soma = NeuritesSoma(0,0,[cx,cy]);
    tree = loadCSVTree(cellfile, soma, scale);
    samplecells{2,m} = tree;
    plotdigineuron(tree,colorstring(m));
    
end
%Find overlapping regions
synapses = findsynapses(samplecells); %2 cells    
sprintf('total synapses=%d', length(synapses))


%%%%METHODS FOR NEURON TREE%%%%%%
%%TODO: DIGITAL NEURON: loop through other tree until find intersect C = intersect(A,B,'rows')
%[xi,yi] = polyxpoly(x1,y1,x2,y2) returns the intersection points of two polylines in a planar, Cartesian system
function synapses = findsynapses(samples)
    synapses = {};
    ctr = 0;
    tree1 = samples{2,1};
    tree2 = samples{2,2};
    for j=1:length(tree2)
        branches2 = tree2{j,1};
        data2 = loadById(branches2);
        N2 = branches2{1,1}(1,1); %root node
        out2 = [];
        out2 = inOrder(data2, N2, out2);
        %figure
        hold on
        for h=1:length(out2)
            N0 = out2(h);
            x1 = N0.points(:,1);
            y1 = N0.points(:,2);
            sprintf('Cell2 node id=%s', N0.id);
            plot(x1,y1,'r-')
            hold on
            for k=1:length(tree1)
                branches = tree1{k,1};
                data = loadById(branches);
                N = branches{1,1}(1,1); %root node
                out = [];
                out = postOrder(data, N, out);
                arrayfun(@(n) plot(n.points(:,1),n.points(:,2),'g-'),out)
                regions = arrayfun(@(n) length(polyxpoly(x1,y1,n.points(:,1),n.points(:,2))>1),out);
                %if positive match - save N - create CSVSynapse
                %[xi,yi]=polyxpoly(x1,y1,s.cell2.parentnode.points(:,1),s.cell2.parentnode.points(:,2))
                if (~isempty(find(regions>=1)))
                    disp("regions found")
                    [idx] = find(regions>=1);
                    for w=1:length(idx)
                        
                        p = out(idx(w));
                        [xi,yi]=polyxpoly(x1,y1,p.points(:,1),p.points(:,2));
                        %TODO split for multiple synapses - each separate
                        for f=1:length(xi)
                            ctr =ctr+1;
                            synapses{ctr}=struct('x',xi(f),'y',yi(f),'cell1', p, 'cell2', N0);
                        end
                        %figure
                        %plot(x1,y1,'g-')
                        
                        %plot(p.points(:,1),p.points(:,2),'r-')
                        plot(xi,yi,'xm','LineWidth',2)
                        
                    end
                end
            end
            %hold off
        end
    end
    
end

%generate missing soma and calculate centroid
function [cx,cy,area,perim,x,y] = findsomacentroid(tree1)
    polycoordsX = [];
    polycoordsY = [];
    for k=1:length(tree1)
        N = tree1{k,1}{1,1}; %root node
        if (N.nodelevel == 1)
            polycoordsX(end+1) =N.points(1,1);
            polycoordsY(end+1) =N.points(1,2);
            
        end
    end
   
    %Sequence with reference to centroid
    x = polycoordsX;
    y = polycoordsY;
    xm = x - mean(polycoordsX);
    ym = y - mean(polycoordsY);
    [t1,r1] = cart2pol(xm,ym)
    mapObj = containers.Map(t1,r1);
    [xp,yp] = pol2cart(cell2mat(mapObj.keys()),cell2mat(mapObj.values()))
    x = xp + mean(polycoordsX);
    y = yp + mean(polycoordsY);
    %Append first set of point to complete polygon
    x(end+1) = x(1);
    y(end+1) = y(1);
    %show in plot
    plot(x,y,'LineStyle',':'); 
    [ geom, iner, cpmo ] = polygeom( x, y );
    area = geom(1);
    cx = geom(2);
    cy = geom(3);
    perim = geom(4);
    %cx = mean(x);
    %cy = mean(y);
    hold on;
    plot(cx,cy,'marker','o','color','b');
end


function plotdigineuron(tree1,color)
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
        colorstring = 'kbgrymc';%colorstring(out(i).nodelevel)
        b = mod(length(colorstring),k)
        if(color=="")
            color = colorstring(b);
        end
        for i=1:length(out)
            plot(out(i).points(:,1),out(i).points(:,2),'color',color); 
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


    