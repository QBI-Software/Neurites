%%Test NeuritesAnalyser and NeuritesSynapse
%GENERATE IMAGE WITH CSV data ONLY
cellfile='sampledata/16-09-15_pair_04_AC.csv'; %010415_DSdata.csv';
colorstring = 'rgkbgrymckbgrymc';%colorstring(out(i).nodelevel)

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
plotdigineuron(tree,'r');
or = 1000/2; %radius
ir = 900/2;
ho = circle(cx,cy,ir,'c');
hi = circle(cx,cy,or,'c');
%Find intersecting regions
regions = findannulusregions(tree,ho,hi); %tree and plots of inner and outer annulus    
sprintf('total regions=%d', length(regions))

function h = circle(x,y,r, color)
    hold on
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(xunit, yunit, color);
    hold off
end

%[xi,yi] = polyxpoly(x1,y1,x2,y2) returns the intersection points of two polylines in a planar, Cartesian system
function regions = findannulusregions(tree, ho,hi)
    regions = {};
    ctr = 0;
   
    for j=1:length(tree)
        branches = tree{j,1};
        data = loadById(branches);
        N = branches{1,1}(1,1); %root node
        out = [];
        out = inOrder(data, N, out);
        %figure
        hold on
        for h=1:length(out)
            N0 = out(h);
            x1 = N0.points(:,1);
            y1 = N0.points(:,2);
            sprintf('Cell node id=%s', N0.id)
            plot(x1,y1,'g-')
            %[ix1,iy1] =polyxpoly(x1,y1,ho.XData,ho.YData);
            %[ix2,iy2] =polyxpoly(x1,y1,hi.XData,hi.YData);
            [ix,iy] =polyxpoly(x1,y1,[ho.XData,hi.XData],[ho.YData,hi.YData]);
            if (~isempty(ix))
                %have intersection with outer and inner
                plot(ix,iy,'mx');
                [in,on] = inpolygon(x1,y1,ix,iy);
                plot(x1(in),y1(in),'y+');
                segx = cat(1,ix(1),x1(in),ix(2));
                segy = cat(1,iy(1),y1(in),iy(2));
                ln = getdistance(segx,segy);
                
                r = sqrt((N0.vol/(pi * N0.branchlength)))
                v = pi * r*r * ln;
                sa = 2 * pi * r * (r + ln);
                ctr =ctr+1;
                regions{ctr}=struct('x',segx,'y',segy,'neurite', N0, 'length',ln,'vol', v, 'sa', sa);
            
            end
            %hold off
        end
    end
    
end

function ln = getdistance(x,y)
    ln = 0;
    for i=1:length(x) - 1 
        ln = ln + pdist([x(i:i+1) y(i:i+1)]);
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


    