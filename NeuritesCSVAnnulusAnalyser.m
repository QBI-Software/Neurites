classdef NeuritesCSVAnnulusAnalyser
    properties
        csvfile     %CSV filename
        tree        %CSV Tree with NeuritesNode objects
        soma        %NeuritesSoma object
        colorstring = 'bmcykbgymc';
        scale = 1;  %allows for adjustment
        regions     %list of matching region x,y coords, NeuritesNode, length, sa, vol
        shape       %[inner diam, outer diam] for annulus, [x,y,width,height] for rectangle
        annulus     %boolean - if annulus or rectangle
        direction   %repeated rectangles in one direction: up,down,right,left
        ctr         %number of matching segments in rois
    end
    methods
        function obj = NeuritesCSVAnnulusAnalyser(csvfile, annulus, shape, scale, direction)
            obj.scale = 1; %scale;
            obj.csvfile = csvfile;
            obj.shape = shape;
            obj.annulus = annulus;
            obj.direction = direction;
            obj.regions ={};
            % Load data via digineuron
            T = readtable(csvfile);
            centroidxy = [mean(T.StartX) mean(T.StartY)]; %estimate from CSV data
            soma = NeuritesSoma(0,0,centroidxy);
            tree = loadCSVTree(csvfile, soma, obj.scale);
            sprintf('Loaded neuron from %s with %d trees', csvfile,length(tree));
            
            %Actual centroid data - reload tree
            [cx,cy,area,perim,x,y] = findsomacentroid(tree);
            obj.soma = NeuritesSoma(area,perim,[cx,cy]);
            obj.tree = loadCSVTree(csvfile, soma, obj.scale);

            plotdigineuron(obj.tree,'r');
            hold on
            obj.ctr = 0;
            if (annulus) % assumes shape with [innerdiam,outerdiam]
                hi = circle(cx,cy,shape(1),'c');
                ho = circle(cx,cy,shape(2),'c');
                [obj.ctr,obj.regions{1}] = obj.findannulusregions(hi,ho);
            else % assumes shape with [x,y,width,height]
                upheight = round(abs(max(T.StartY)));
                downheight = round(abs(min(T.StartY)));
                leftwidth = round(abs(min(T.StartX)));
                rightwidth = round(abs(max(T.StartX)));
                x = cx;
                y = cy;
                width = shape(1);
                height = shape(2);
                
                switch direction
                    case 'up'
                        numrects = ceil(upheight/height)
                        x = x - width/2;
                        for i = 1:numrects 
                            y = y + height;
                            plot(x,y,'mx');
                            xdata = [x,x+width,x+width,x,x];
                            ydata = [y,y,y-height,y-height,y];
                            plot(xdata,ydata,obj.colorstring(i));
                            [ctr,obj.regions{i}] = obj.findrectangleregions(xdata,ydata);
                            obj.ctr = obj.ctr + ctr;
                        end

                    case 'down'
                        numrects = ceil(downheight/height)
                        x = x - width/2;
                        for i = 1:numrects 
                            plot(x,y,'mx');
                            xdata = [x,x+width,x+width,x,x];
                            ydata = [y,y,y-height,y-height,y];
                            plot(xdata,ydata,obj.colorstring(i));
                            [ctr,obj.regions{i}] = obj.findrectangleregions(xdata,ydata);
                            obj.ctr = obj.ctr + ctr;
                            y = y - height;
                        end
                    case 'left'
                        numrects = ceil(leftwidth/width)
                        y = y + height/2;
                        for i = 1:numrects 
                            x = x - width;
                            plot(x,y,'mx');
                            xdata = [x,x+width,x+width,x,x];
                            ydata = [y,y,y-height,y-height,y];
                            plot(xdata,ydata,obj.colorstring(i));
                            [ctr,obj.regions{i}] = obj.findrectangleregions(xdata,ydata);
                            obj.ctr = obj.ctr + ctr;
                            
                        end
                    case 'right'
                        numrects = ceil(rightwidth/width)
                        y = y + height/2;
                        for i = 1:numrects 
                            
                            plot(x,y,'mx');
                            xdata = [x,x+width,x+width,x,x];
                            ydata = [y,y,y-height,y-height,y];
                            plot(xdata,ydata,obj.colorstring(i));
                            [ctr,obj.regions{i}] = obj.findrectangleregions(xdata,ydata);
                            obj.ctr = obj.ctr + ctr;
                            x = x + width;
                        end
                    case 0 %all dims provided
                        if (length(shape) == 4)
                            x = shape(1);
                            y = shape(2);
                            width = shape(3);
                            height = shape(4);
                            plot(x,y,'mx');
                            xdata = [x,x+width,x+width,x,x];
                            ydata = [y,y,y-height,y-height,y];
                            plot(xdata,ydata,'m');
                            [ctr,obj.regions{1}] = obj.findrectangleregions(xdata,ydata);
                            obj.ctr = obj.ctr + ctr;
                        else %centred around soma
                            x = cx - width/2;
                            y = cy + height/2; 
                            plot(x,y,'mx');
                            xdata = [x,x+width,x+width,x,x];
                            ydata = [y,y,y-height,y-height,y];
                            plot(xdata,ydata,'m');
                            [ctr,obj.regions{1}] = obj.findrectangleregions(xdata,ydata);
                            obj.ctr = obj.ctr + ctr;
                        end

                end
                 
            end
        end

       

        %[xi,yi] = polyxpoly(x1,y1,x2,y2) returns the intersection points of two polylines in a planar, Cartesian system
        function [ctr,regions] = findannulusregions(obj, hi,ho)
            
            regions = {};
            ctr = 0;

            for j=1:length(obj.tree)
                branches = obj.tree{j,1};
                data = loadById(branches);
                N = branches{1,1}(1,1); %root node
                out = [];
                out = inOrder(data, N, out);
                %filter those within range soma dist
                filteredout = arrayfun(@(n) compareDistanceToSoma(obj.soma.centroid, n.points,obj.shape), out);
                out = out(filteredout);
                %figure
                hold on
                for h=1:length(out)
                    N0 = out(h);
                    x1 = N0.points(:,1);
                    y1 = N0.points(:,2);
                    sprintf('Cell node id=%d, h=%d', N0.id,h)
                    %plot(x1,y1,'g-');
                    %[ix,iy] =polyxpoly(x1,y1,[ho.XData,hi.XData],[ho.YData,hi.YData]);
                    [ix1,iy1,ii1] =polyxpoly(x1,y1,hi.XData,hi.YData); %inner match
                    [ix2,iy2,ii2] =polyxpoly(x1,y1,ho.XData,ho.YData); %outer match
                    ix = [ix1,ix2];
                    iy = [iy1,iy2];
                    multiseg ={};
                    if (~isempty(ix))
                        if(length(ix) == 2)
                            %have intersection with outer and inner
                            plot(ix,iy,'mx');
                            %[in,on] = inpolygon(x1,y1,ix,iy); %intermittent fails
                            [in1,on1] = inpolygon(x1,y1,ho.XData,ho.YData); % all in outer circle
                            [in2,on2] = inpolygon(x1,y1,hi.XData,hi.YData); % all in inner circle 
                            in = xor(in1,in2);
                            %plot(x1(in),y1(in),'go');
                            segx = cat(1,ix(1),x1(in),ix(2));
                            segy = cat(1,iy(1),y1(in),iy(2));
                            multiseg{1} = struct('x',segx,'y',segy);
                        elseif(length(ix) > 2)
                            if (~isempty(ii1))
                                multiseg = findsubsegments(multiseg,ix1,iy1,x1,y1,hi.XData,hi.YData);
                            else
                                multiseg = findsubsegments(multiseg,ix2,iy2,x1,y1,ho.XData,ho.YData);
                            end
                            
                        else %only one intersection 
                            %direction detect
                            fwd = getDirectionNeurite(obj.soma.centroid,x1,y1);
                            p = pdist2([ix,iy],[x1,y1]);
                            idx = find(p == min(p));
                            
                            %detect inner or outer
                            %[ix1,iy1] =polyxpoly(x1,y1,hi.XData,hi.YData);
                            
                            if ((~isempty(ix1) && fwd) || (isempty(ix1) && ~fwd))  
                                    segx = cat(1, ix(1), x1(idx+1:end));
                                    segy = cat(1, iy(1), y1(idx+1:end));
                            else %add to end
                                    segx = cat(1, x1(1:idx), ix(1));
                                    segy = cat(1, y1(1:idx), iy(1));
                            end
                            multiseg{1} = struct('x',segx,'y',segy);        
                         
                        end
                    else % no intersection
                        segx = x1;
                        segy = y1;
                        multiseg{1} = struct('x',segx,'y',segy);
                    end         
                    
                    for i=1:length(multiseg)
                        segx = multiseg{i}.x;
                        segy = multiseg{i}.y;
                        plot(segx,segy,'m+');
                        ln = getdistance(segx,segy);
                        r = sqrt((N0.vol/(pi * N0.branchlength)));
                        v = pi * r*r * ln;
                        sa = 2 * pi * r * (r + ln);
                        ctr =ctr+1;
                        regions{ctr}=struct('x',segx,'y',segy,'neurite', N0, 'length',ln,'vol', v, 'sa', sa);
                    end
                   
                    %hold off
                end
            end

        end
        
        function [ctr,regions] = findrectangleregions(obj,xdata,ydata)
            regions = {};
            ctr = 0; %length(obj.regions);

            for j=1:length(obj.tree)
                branches = obj.tree{j,1};
                data = loadById(branches);
                N = branches{1,1}(1,1); %root node
                out = [];
                out = inOrder(data, N, out);
                hold on
                for h=1:length(out)
                    N0 = out(h);
                    x1 = N0.points(:,1);
                    y1 = N0.points(:,2);
                    sprintf('Cell node id=%d, h=%d', N0.id,h)
                    plot(x1,y1,'g-');
                    [ix,iy] = polyxpoly(x1,y1,xdata,ydata); % edges
                    [in,on] = inpolygon(x1,y1,xdata,ydata); % contents
                    multiseg ={};

                    if (~isempty(in(in==true)))
                        switch length(ix)
                            case 2
                                segx = cat(1,ix(1),x1(in),ix(2));
                                segy = cat(1,iy(1),y1(in),iy(2));
                                multiseg{1} = struct('x',segx,'y',segy);
                            case 1
                                %prepend or postpend
                                fwd = getDirectionNeurite(obj.soma.centroid,x1,y1);
                                p = pdist2([ix,iy],[x1(in),y1(in)]);
                                idx = find(p == min(p));

                                if ((idx <= 1 && fwd) || (idx >= length(x1(in)) && ~fwd))  
                                        segx = cat(1, ix(1), x1(in));
                                        segy = cat(1, iy(1), y1(in));
                                else %add to end
                                        segx = cat(1, x1(in), ix(1));
                                        segy = cat(1, y1(in), iy(1));
                                end
                                multiseg{1} = struct('x',segx,'y',segy);
                            case 0
                                segx = x1(in);
                                segy = y1(in);
                                multiseg{1} = struct('x',segx,'y',segy);
                            otherwise
                                %multiseg = findsubsegments(multiseg,ix,iy,x1,y1,xdata,ydata);
                                p = pdist2([ix,iy],[x1,y1]);
                                s = size(p);
                                mc = length(multiseg) + 1;
                                start = 0;
                                segs =[];
                                for i=1:s(1)
                                    r = p(i,:); %per row
                                    idx = find(r == min(r)) %find idx of least distance from neurite
                                    segs(i) = idx;
                                end
                                sortsegs = sort(unique(segs));
                                for i=1:length(sortsegs)
                                    %Split neurite into segs
                                    segx = x1(start+1:sortsegs(i)-1);
                                    segy = y1(start+1:sortsegs(i)-1);
                                    [in,on] = inpolygon(segx,segy,xdata,ydata);
                                    if (~isempty(in(in==true)))
                                        segx = x1(start:sortsegs(i));
                                        segy = y1(start:sortsegs(i));
                                        multiseg{mc} = struct('x',segx,'y',segy);
                                    end
                                    start = sortsegs(i);
                                end         
                        end
                       
                        for i=1:length(multiseg)
                            segx = multiseg{i}.x;
                            segy = multiseg{i}.y;
                            plot(segx,segy,'m+');
                            ln = getdistance(segx,segy);
                            r = sqrt((N0.vol/(pi * N0.branchlength)))
                            v = pi * r*r * ln;
                            sa = 2 * pi * r * (r + ln);
                            cx = xdata(1) + abs((xdata(2)- xdata(1))/2);
                            cy = ydata(1) - abs((ydata(3) - ydata(1))/2);
                            ctr =ctr+1;
                            regions{ctr}=struct('x',segx,'y',segy,'neurite', N0, 'length',ln,'vol', v, 'sa', sa, 'cx',cx,'cy',cy);
                        end
                    end
                end
            end
        end

        
        %% Generate Table data - export data to a table format
        function [colnames,tabledata] = generateTable(obj)
            colnames = {'Region' 'Cx' 'Cy' 'Area' ...
                'Tree' 'Order' 'Branchlength' 'Rownum'  ...
                'Branch_angle' 'Branch_vol' 'Branch_sa' ...
                'Seg_length' 'Seg_vol' 'Seg_sa' };
           
            c = 1;
            tabledata =zeros(length(obj.ctr),length(colnames)); %preallocate
            if (obj.annulus) 
                area = ((obj.shape(2))^2 * pi) - ((obj.shape(1))^2 * pi);
                for i=1:length(obj.regions)
                    region = obj.regions{1,i};
                    for j=1:length(region)
                        r = region{1,j};
                        row1 = [i obj.soma.centroid(1) obj.soma.centroid(2) area ...
                            r.neurite.tree r.neurite.nodelevel r.neurite.branchlength r.neurite.id ...
                            r.neurite.anglearc r.neurite.vol r.neurite.sa ...
                            r.length r.vol r.sa];
                        tabledata(c,:) = row1;
                        c = c+1;
                    end

                end
            else %rectangles
                
                if (length(obj.shape) ==2)
                    area = obj.shape(1) * obj.shape(2);
                else
                    area = obj.shape(3) * obj.shape(4);
                end
               
                for i=1:length(obj.regions)
                    region = obj.regions{1,i};
                    for j=1:length(region)
                        r = region{1,j};
                        row1 = [i r.cx r.cy area ...
                            r.neurite.tree r.neurite.nodelevel r.neurite.branchlength r.neurite.id ...
                            r.neurite.anglearc r.neurite.vol r.neurite.sa ...
                            r.length r.vol r.sa];
                        tabledata(c,:) = row1;
                        c = c+1;
                    end

                end
            end
            tabledata( ~any(tabledata,2), : ) = [];  %remove zero rows
        end
        
        %Output segment data corresponding to original CSV file
        function [T1,T2] = outputSegmentXY(obj)
            T1 = readtable(obj.csvfile); 
            %load row idxs
            idxs = [];
            %collect all xy pts - may be duplicates
            allpts = [];
            regions ={};
           
            
            for i=1:length(obj.regions)
                region = obj.regions{1,i};
                    for j=1:length(region)
                        r = region{1,j};
                        %idxs = cat(1, idxs,r.neurite.id); %this is only branch last pt
                        %fR = find((StartY >= yC-tol) & (StartY <=yC+tol));
                        a = find(ismember(r.neurite.points,r.x));
                        ais = r.neurite.id - (length(r.neurite.points) - a);
                        if (ais > 0)
                            idxs = cat(1, idxs,ais);
                            allpts = cat(1, allpts,[r.x r.y]);
                        end
                    end

            end
            T1 = T1(idxs,:);
            T2 = table(allpts(:,1),allpts(:,2),'VariableNames',{'X' 'Y'});
        end
        
     
    end %end class methods
end %end class

%calculate distance point by point assuming straight lines between points
function ln = getdistance(x,y)
    ln = 0;
    for i=1:length(x) - 1 
        ln = ln + pdist([x(i:i+1) y(i:i+1)]);
    end
end

%determine if neurite is radiating outward or wrapping back towards soma
function fwd = getDirectionNeurite(centroid, x1,y1)
    s0 = pdist2(centroid,[x1(1),y1(1)]);
    s1 = pdist2(centroid,[x1(end),y1(end)]);
    %forward
    fwd = (s0 < s1);
end

 function h = circle(x,y,r, color)
   % hold on
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(xunit, yunit, color);
    %hold off
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
if (length(polycoordsX) > 2)
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
else
    cx = mean(polycoordsX);
    cy = mean(polycoordsY);
    area = 0;
    perim = 0;
end

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

%compare distance to soma of begin and end of points greater than shape
% [innerdiam,outerdiam]
function val = compareDistanceToSoma(centroid,points, shape)
    id = shape(1); %innerdiam
    od = shape(2); %outerdiam
    p1 = pdist2(centroid,points(1,:)); %distance from start of points
    p2 = pdist2(centroid,points(end,:)); %distance from end of points
    val1 = ((p1 > id) || (p2 > id));
    val2 = (p1 > od && p2 > od);
    val = xor(val1,val2); %true if either are true but not both
end

function multiseg = findsubsegments(multiseg,ix,iy,x1,y1,xdata,ydata)
    p =pdist2([ix,iy],[x1,y1]);
    s = size(p);
    mc = length(multiseg) + 1;
    cache = 0;
    for i=1:s(1)
        r = p(i,:); %per row
        idx = find(r == min(r)) %find idx of least distance from line
        if (~cache)
            x1c = x1(idx+1:end);
            y1c = y1(idx+1:end); 
            segx = cat(1, x1(idx), x1c);
            segy = cat(1, y1(idx), y1c);
        else
            x1c = x1(idx+1:cache-1);
            y1c = y1(idx+1:cache-1);
            segx = cat(1, x1(idx), x1c, x1(cache));
            segy = cat(1, y1(idx), y1c, y1(cache));
        end
        [in,on] = inpolygon(x1c,y1c,xdata,ydata); %test segment without matching pts

        if (length(in(in==true))> 0)
            plot(segx,segy,'y','LineWidth',2);
            multiseg{mc} = struct('x',segx,'y',segy);
            mc = mc + 1;
        end
        cache = idx    
    end
    % check final seg
    x1c = x1(1:idx-1);
    y1c = y1(1:idx-1);
    [in,on] = inpolygon(x1c,y1c,xdata,ydata);
    if (length(in(in==true))> 0)
        segx = cat(1, x1c, x1(idx));
        segy = cat(1, y1c, y1(idx));
        %plot(segx,segy,'y','LineWidth',2);
        multiseg{mc} = struct('x',segx,'y',segy);
    end
end




