classdef NeuritesCSVAnalyser
    properties
        samplecells %cellarray of csv filenames(samplecells{1,:}), trees (samplecells{2,:}) and somas(samplecells{3,:})
        colorstring = 'rgkbgrymckbgrymc';
        scale = 1; %default if using csv data
        synapses %list of synapse x,y coords, with cell1 and cell2 as NeuritesNodes where points are crossing
    end
    methods
        function obj = NeuritesCSVAnalyser(csv1file,csv2file, scale)
            obj.scale = scale;
            obj.samplecells = {csv1file, csv2file};
            % Load data via digineuron
            for m=1:length(obj.samplecells)
                cellfile = obj.samplecells{1,m};
                T = readtable(cellfile);
                centroidxy = [mean(T.StartX) mean(T.StartY)]; %estimate from CSV data
                %scale = 1;
                soma = NeuritesSoma(0,0,centroidxy);
                tree = loadCSVTree(cellfile, soma, obj.scale);
                sprintf('Loaded neuron from %s with %d trees', cellfile,length(tree));
                %Actual centroid data - reload tree
                [cx,cy,area,perim,x,y] = findsomacentroid(tree);
                soma = NeuritesSoma(area,perim,[cx,cy]);
                tree = loadCSVTree(cellfile, soma, obj.scale);
                obj.samplecells{2,m} = tree;
                obj.samplecells{3,m} = soma;
                plotdigineuron(tree,obj.colorstring(m));
                hold on

            end
            obj = obj.findsynapses();
        end


        function obj = set.synapses(obj,segments)
            if (~isempty(segments))
                obj.synapses = segments;
            else
                error('Error: Cannot set SYNAPSES - Empty set')
            end
        end

        function obj = findsynapses(obj)
            list_synapses = {}; %init
            ctr = 0;
            tree1 = obj.samplecells{2,1};
            tree2 = obj.samplecells{2,2};
            for j=1:length(tree2)
                branches2 = tree2{j,1};
                data2 = loadById(branches2);
                N2 = branches2{1,1}(1,1); %root node
                out2 = [];
                out2 = inOrder(data2, N2, out2);
                %figure
                %hold on
                for h=1:length(out2)
                    N0 = out2(h);
                    x1 = N0.points(:,1);
                    y1 = N0.points(:,2);
                    %sprintf('Cell2 node id=%s', N0.id);
                    %plot(x1,y1,'r-')
                    %hold on
                    for k=1:length(tree1)
                        branches = tree1{k,1};
                        data = loadById(branches);
                        N = branches{1,1}(1,1); %root node
                        out = [];
                        out = postOrder(data, N, out);
                        %arrayfun(@(n) plot(n.points(:,1),n.points(:,2),'g-'),out)
                        regions = arrayfun(@(n) length(polyxpoly(x1,y1,n.points(:,1),n.points(:,2))>1),out);
                        %if positive match - save N - create CSVSynapse
                        %[xi,yi]=polyxpoly(x1,y1,s.cell2.parentnode.points(:,1),s.cell2.parentnode.points(:,2))
                        if (~isempty(find(regions>=1)))
                            %disp("regions found")
                            [idx] = find(regions>=1);
                            for w=1:length(idx)

                                p = out(idx(w));
                                [xi,yi]=polyxpoly(x1,y1,p.points(:,1),p.points(:,2));
                                %Split for multiple synapses - each separate
                                for f=1:length(xi)
                                    ctr =ctr+1;
                                    list_synapses{ctr}=struct('x',xi(f),'y',yi(f),'cell1', p, 'cell2', N0);
                                end
                                plot(xi,yi,'xm','LineWidth',2)

                            end
                        end
                    end
                    %hold off
                end
            end
            %Save data
            obj.synapses = list_synapses;
        end

        function d = getDistanceToSoma(obj,synapsenode,x,y)
            %subtract synapse from end of segment length
            synxy = [x,y];
            endxy = [synapsenode.points(end,1), synapsenode.points(end,2)];
            d = -abs((distanceBetween(synxy,endxy,'euclidean')));

            while(~isempty(synapsenode.parentnode))
                d = d+ synapsenode.branchlength;
                synapsenode = synapsenode.parentnode;
            end

        end
        %% Generate Table data - export data to a table format
        function [colnames,tabledata] = generateTable(obj,cell1label, cell2label)
            colnames = {'Synapse_X' 'Synapse_Y' ...
                'Cell1_tree' 'Cell1_order' 'Cell1_branchlength' 'Cell1_distance2soma'  ...
                'Cell1_theta' 'Cell1_rho' 'Cell1_deg' ...
                'Cell2_tree' 'Cell2_order' 'Cell2_branchlength' 'Cell2_distance2soma'  ...
                'Cell2_theta' 'Cell2_rho' 'Cell2_deg' };
            colnames = strrep(colnames, 'Cell1', cell1label);
            colnames = strrep(colnames, 'Cell2', cell2label);

            tabledata =zeros(length(obj.synapses),length(colnames)); %preallocate
            soma1 = obj.samplecells{3,1};
            soma2 = obj.samplecells{3,2};

            for i=1:length(obj.synapses)
                syn = obj.synapses{1,i};
                d1 = obj.getDistanceToSoma(syn.cell1,syn.x,syn.y);
                d2 = obj.getDistanceToSoma(syn.cell2,syn.x,syn.y);
                [theta1,rho1,deg1] = findAngleSoma(syn.x,syn.y,soma1,obj.scale);
                [theta2,rho2,deg2] = findAngleSoma(syn.x,syn.y,soma2,obj.scale);
                row1 = [syn.x syn.y ...
                    syn.cell1.tree syn.cell1.nodelevel syn.cell1.branchlength d1 ...
                    theta1 rho1 deg1 ...
                    syn.cell2.tree syn.cell2.nodelevel syn.cell2.branchlength d2 ...
                    theta2 rho2 deg2 ];

                %tabledata = cat(1,tabledata,row1);
                tabledata(i,:) = row1;

            end
            tabledata( ~any(tabledata,2), : ) = [];  %remove zero rows
        end

    end %end class methods
end %end class

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
[t1,r1] = cart2pol(xm,ym);
mapObj = containers.Map(t1,r1);
[xp,yp] = pol2cart(cell2mat(mapObj.keys()),cell2mat(mapObj.values()));
x = xp + mean(polycoordsX);
y = yp + mean(polycoordsY);
%Append first set of point to complete polygon
x(end+1) = x(1);
y(end+1) = y(1);

%calculate geometric params
[ geom, iner, cpmo ] = polygeom( x, y );
area = geom(1);
cx = geom(2);
cy = geom(3);
perim = geom(4);
%show in plot
plot(x,y,'LineStyle',':');
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
    colorstring = 'grbymc';%colorstring(out(i).nodelevel)
    b = mod(length(colorstring),k) + 1;
    if(isempty(color))
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




