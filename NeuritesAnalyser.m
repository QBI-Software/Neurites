classdef NeuritesAnalyser
   properties
      I
      Iroi
      rbw
      gbw
      maskedR
      maskedG
      Segments
      Synapses
      
   end
   methods
      function obj = NeuritesAnalyser(rgbimg,roi,R,G)
         if (ischar(rgbimg))
          obj.I = imread(rgbimg);
         else
          obj.I = rgbimg;
         end
         if (ischar(roi))
          obj.Iroi = imread(roi);
         else
          obj.Iroi = roi;
         end
         %channel separation if split cells not supplied as bw
         if (~exist('R', 'var') || isempty(R))
            R = obj.I(:,:,2);
            obj.rbw = im2bw(R);
            obj.rbw = imcomplement(obj.rbw);
         else
            obj.rbw = R;
         end
         if (~exist('G', 'var') ||isempty(G))
             G = obj.I(:,:,1);
            obj.gbw = im2bw(G);
            obj.gbw = imcomplement(obj.gbw);
         else
            obj.gbw = G;
         end
         

         obj.maskedR = obj.rbw; % Simply a copy at first.
         %attempt to reduce to skeleton
         %obj.maskedR = bwmorph(obj.rbw,'skel', Inf);
         obj.maskedR(~obj.Iroi) = 0;  % Set all non-keeper pixels to zero.

         obj.maskedG = obj.gbw; % Simply a copy at first.
         %obj.maskedG = bwmorph(obj.gbw,'skel', Inf);
         obj.maskedG(~obj.Iroi) = 0;  % Set all non-keeper pixels to zero.
         
      end
      
      function obj = set.Segments(obj,segments)
          if (~isempty(segments))
             obj.Segments = segments;
          else
             error('Error: Cannot set segments - Empty set')
          end
      end
      function obj = set.Synapses(obj,segments)
          if (~isempty(segments))
             obj.Synapses = segments;
          else
             error('Error: Cannot set synapses - Empty set')
          end
      end
      %Processes image to find binary images: 
      % skeleton (BWS), branchpoints (BWB), endpoints (BWE)
      % and list of lines (lines)
      function [lines,bp,ep] = processMorph(obj,rbw)
       
        %Smoothing elements
        se = strel('line',2,5);
        se1 = strel('line',3,0);
        se2 = strel('line',3,90);
        %rbw = bwmorph('clean');
        %BW3 = bwmorph(rbw,'open');
        BW1 = imdilate(rbw,[se1 se2],'full');
        BW2 = bwmorph(BW1,'bridge');
        BW3 = imerode(BW2, se);
        BWS = bwmorph(BW3,'skel',Inf);
        BWB = bwmorph(BWS,'branchpoints');
        BWE = bwmorph(BWS,'endpoints');
        % get rid of overlap between ep and bp
        BWE2 = bwmorph(BWE, 'thicken');
        BWB2 = bwmorph(BWB, 'thicken');
        BWB(find(BWE2)) = 0;
        BWE(find(BWB2)) = 0; 
        %Plots for checking
        figure, imshow(rbw); hold on;
        [m,n] = find(BWB == 1);
        %reduce by distance
        D = pdist2(m,n,'euclidean');
        [r,c]= find(D <= 3);
        %removes duplicates within range but should average
        BWB(ismember(r,c),:)=0;
        [m,n] = find(BWB == 1);
        bp = {length(n)};
        for i=1:length(n)
            plot(n(i),m(i),...
                 'color','b','marker','X','linestyle','none','LineWidth',2);
            bp{i} = [n(i),m(i)];
        end
        
        [j,k] = find(BWE == 1);
        lines = {length(k)};
        ctr = 1;
        ep = {length(k)};
        for i=1:length(k)
            r = j(i);
            c = k(i);
            contour = bwtraceboundary(BWS,[r c],'S',8,400,'counterclockwise');
            if (~isempty(contour))
                plot(contour(:,2),contour(:,1),'g','LineWidth',2);
                lines{ctr} = contour;
                ctr = ctr+1;
            end
            plot(k(i),j(i),...
                 'color','r','marker','O','linestyle','none','LineWidth',2);
            ep{i} = [n(i),m(i)];
        end
        legend('x - Branchpoints','o - Endpoints')
        hold off
        
      end
      
      function obj = findSegments(obj,minlength,tol)
          Rb = bwboundaries(obj.maskedR);
          Gb = bwboundaries(obj.maskedG);
          %[Rb,Rbp,Rep] = obj.processMorph(obj.rbw);
          %[Gb,Gbp,Gep] = obj.processMorph(obj.gbw);
          
          nR = length(Rb)
          nG = length(Gb)
          match = 1;
          segments = {};
          resolution = 1;
          minsynapselength=4;
          %test
          %figure
          imshow(obj.maskedR)
          hold on
          %Masked region coords
          [mar,mbr] = find(obj.maskedR==1);
          [mag,mbg] = find(obj.maskedG==1);
          for i = 1:nR
            redSegment = Rb{i};
            if (isempty(find(ismember(redSegment,[mar,mbr])==1, 1)))
                plot(redSegment(:,2),redSegment(:,1),'b','LineWidth',2)
                continue
            else
                plot(redSegment(:,2),redSegment(:,1),'r','LineWidth',2)
            end
            if (length(redSegment) > minlength)
                x1 = redSegment(:,2);
                y1 = redSegment(:,1);
                for k = 1:nG
                    greenSegment = Gb{k};
                    if (isempty(find(ismember(greenSegment,[mag,mbg])==1,1)))
                       % plot(greenSegment(:,2),greenSegment(:,1),'y','LineWidth',2)
                        continue
                    end
                    x2 = greenSegment(:,2);
                    y2 = greenSegment(:,1);
                    nx1 = size(x1,1);
                    nx2 = size(x2,1);

                    % calculate all pairwise distances
                    %find gives indices which match values
                    K = find(x2 >= (min(x1) - tol) & x2 <= (max(x1) + tol));
                    O = find(y2 >= (min(y1) - tol) & y2 <= (max(y1) + tol));
                    %find matching indices in x2,y2
                    if (~isempty(K) && ~isempty(O))
                        changes = zeros(nx1,nx2);
                        for m = 1:nx2
                            %changes(:,m) = sqrt((x1-x2(m)).^2 + (y1-y2(m)).^2);
                            changes(:,m) = (abs(x1-x2(m)) + abs(y1-y2(m)))/2;
                        end
                        [j,n] = find(changes <= tol);
                        if (~isempty(j))

                            %check that range is single line
                            blankimg = zeros(size(obj.rbw));
                            blankimg(y1(j),x1(j)) = 1;
                            blankimg(y2(n),x2(n)) = 1;
                            [L, num] = bwlabel(blankimg, 8);

                            if (num) > 1
                                toggle=logical(1);
                                for w=1:num
                                    L1 = L == w;
                                    [col,row] = find(L1);
                                    %L1 = y2(col),x2(row)?
                                    s  = regionprops(L1, 'centroid');
                                    
                                     if toggle
                                        segments{match}{1} = redSegment;
                                        segments{match}{2} = greenSegment;
                                        if (~inpolygon(s.Centroid(1),s.Centroid(2),redSegment(:,2),redSegment(:,1)))
                                            segments{match}{3} = [col,row];
                                            segments{match}{5} = findLength([col,row], resolution);
                                
                                        elseif (~inpolygon(s.Centroid(1),s.Centroid(2),greenSegment(:,2),greenSegment(:,1)))
                                            segments{match}{4} = [col,row];
                                            segments{match}{6} = findLength([col,row], resolution);
                                            toggle = ~toggle;
                                            match = match+1;
                                        end
                                     end
                                end
                            else
                                % CHANGE TO CREATE SYNAPSE
                                segments{match}{1} = redSegment;
                                segments{match}{2} = greenSegment;
                                segments{match}{3} = [y1(j),x1(j)];
                                segments{match}{4} = [y2(n),x2(n)];
                                segments{match}{5} = findLength([y1(j),x1(j)], resolution);
                                segments{match}{6} = findLength([y2(n),x2(n)], resolution);
                                match = match + 1;
                            end
                            
                        end
                        %one match per redsegment
                        continue
                    end
                end
            end
          end
          
          hold off
          if (~isempty(segments))
              obj.Segments = segments;
              status = sprintf('Matched=%d of %d', match, nG)
          else
              status = 'No segments found'
          end
   end
        
    function obj = analyseSynapses(obj,showplots,showtypes)
        
        ctr = 1;
        if (showplots)
            imshow(obj.maskedR)
            hold on;
        end
        synapses = {};
        for p=1:length(obj.Segments)
            type = 1; %en-passant as default
            redSegment = obj.Segments{p}{1};
            greenSegment = obj.Segments{p}{2};
            redhighlight = obj.Segments{p}{3};
            greenhighlight = obj.Segments{p}{4};
            if (length(redhighlight)==0)
                continue
            elseif(length(greenhighlight)==0)
                continue
            end
            
            plot(redSegment(:,2), redSegment(:,1), 'r', 'LineWidth', 2);
            plot(greenSegment(:,2), greenSegment(:,1),'g','LineWidth',2);
            plot(redhighlight(:,2), redhighlight(:,1), 'y', 'LineWidth', 3);
            plot(greenhighlight(:,2), greenhighlight(:,1),'y','LineWidth',3);

            %find min max points on each line
            [p1,p2] = findLinePoints(redhighlight);
            [p3,p4] = findLinePoints(greenhighlight);

            % p4 not used?
            [a12,a13,a23] = findTriangleAngles(p1,p2,p3);

            %ignore if isosceles
            a = [round(a12,0) round(a13,0) round(a23,0)];
            a = sort(a);
            b = [45 45 90];
            c = [0 0 180];
            if (isequaln(a,b) || isequaln(a,c) || ~isempty(a(isnan(a))))
                type = 3; %isosceles indicates intersection not synapse
            elseif(~isempty(intersect(a,b))) %if one angle is present check with other point
                [a12,a13,a23] = findTriangleAngles(p1,p2,p4);
                a = [round(a12,0) round(a13,0) round(a23,0)];
                a = sort(a);
                if (isequaln(a,b))
                    type = 3;
                end
            end
            
            %determine median synapse location
            rmarkerx = median(redhighlight(:,2));
            rmarkery = median(redhighlight(:,1));
            gmarkerx = median(greenhighlight(:,2));
            gmarkery = median(greenhighlight(:,1));
            nrange = [1 1];
            %endpoint of one cell on another
            if (isEndpoint(rmarkerx, rmarkery, redSegment,nrange) || ...
                isEndpoint(gmarkerx,gmarkery,greenSegment,nrange))
                type = 2;
            end
            mycolors = ['m','b','c'];
            if (showplots && ismember(type,showtypes))
                plot(p1(1), p1(2), p2(1), p2(2), p3(1), p3(2),...
                'color',mycolors(type),'marker','o','LineStyle','-','LineWidth', 2);
                plot(rmarkerx,rmarkery,'color',mycolors(type),'marker','X','linestyle','none','LineWidth', 2);
                plot(gmarkerx,gmarkery,'color',mycolors(type),'marker','X','linestyle','none','LineWidth', 2);
            end

            %Save synapses to list
            syn = NeuritesSynapse(redSegment, greenSegment, redhighlight, greenhighlight, ...
               rmarkerx,rmarkery,gmarkerx,gmarkery,type);
            synapses{ctr} = syn;
            ctr = ctr + 1;

        end
        obj.Synapses = synapses;
        if (showplots)
            hold off;
        end
    end
    
    function obj = measureSynapses(obj,showtypes,cell1file,cell2file,scale,shiftx,shifty,tol)
        I = obj.rbw;
        T1 = readtable(cell1file);
        T2 = readtable(cell2file);
        %variable shift TODO adjust
        %shiftx = 177;
        %shifty = -942;
        %plot
        X = (T1.StartX * scale) + shiftx;
        Y = (T1.StartY * scale) + shifty;
        figure
        imshow(I);
        hold on;
        %plot(synR(1),synR(2),'Marker', '*', 'color','m');
        %plot(synG(1),synG(2),'Marker', '*', 'color','y');
        plot(X,-Y,'color','c','LineStyle','-','LineWidth', 2);
        hold off;
        %tol = 10;
        for i = 1: length(obj.Synapses)
            syn = obj.Synapses{i};
            if (ismember(syn.SynapseType,showtypes))
                %reset
                fR = [];
               
                %C1 and C2 median should be similar
                if (isempty(find(abs(syn.MedianC1-syn.MedianC2))> tol))
                    x = mean(syn.MedianC1(1), syn.MedianC2(1));
                    y = mean(syn.MedianC1(2), syn.MedianC2(2));
                else
                    x = syn.MedianC1(1);
                    y = syn.MedianC1(2);
                end
                
                %cell1
                syn.shiftx = shiftx;
                syn.shifty = shifty;
                syn.scale = scale;
                [xC, yC]= syn.img2Coords(x,y);
                %select idx of row containing coords in csv - interpolation
                fR=findCSVIndex(xC,yC,T1.StartX,T1.StartY,tol);
                if (isempty(fR))
                    warning('Warning: C1 synapse coords not found - %d, %d', xC,yC)
                    %continue
                    syn.StartXY1 = [0 0];
                    syn.NeuriteLengthC1 = 0;
                    syn.DistanceC1= 0;
                    syn.BranchPointC1 = 0;
                    syn.BranchTypeC1 = 0;
                else
                    %load it up
                    syn.StartXY1 = [T1.StartX(fR), T1.StartY(fR)];
                    syn.NeuriteLengthC1 = T1.Length__m_(fR)
                    syn.DistanceC1= T1.LengthToBeginning__m_(fR)
                    syn.BranchPointC1 = T1.Order(fR)
                    syn.BranchTypeC1 = T1.PointType(fR)
                end
                %cell2
                fR=findCSVIndex(xC,yC,T2.StartX,T2.StartY,tol);
                if (isempty(fR))
                    warning('Warning: C2 synapse coords not found - %d, %d', xC,yC)
                    %continue
                    syn.StartXY2 = [0 0];
                    syn.NeuriteLengthC2 = 0;
                    syn.DistanceC2= 0;
                    syn.BranchPointC2 = 0;
                    syn.BranchTypeC2 = 0;
                else
                    %load it up
                    syn.StartXY2 = [T2.StartX(fR), T2.StartY(fR)];
                    syn.NeuriteLengthC2 = T2.Length__m_(fR)
                    syn.DistanceC2= T2.LengthToBeginning__m_(fR)
                    syn.BranchPointC2 = T2.Order(fR)
                    syn.BranchTypeC2 = T2.PointType(fR)
                end
                %save it
                obj.Synapses{i} = syn;
                
            end
        end
       


    end
    function [colnames,tabledata] = generateTable(obj,types)
        colnames = {'Type','Cell1 X','Cell1 Y','Cell2 X','Cell2 Y',...
            'Cell1 length','Cell1 distance','Cell1 order', ...
            'Cell2 length','Cell2 distance','Cell2 order', ...
            'Cell1 StartX','Cell1 StartY','Cell2 StartX','Cell2 StartY' };
        tabledata =[];
        for i=1:length(obj.Synapses)
            syn = obj.Synapses{i};
            if (ismember(syn.SynapseType,types))
                
                %show data in table
                row1 = [syn.SynapseType ...
                    syn.MedianC1 syn.MedianC2 ...
                    syn.NeuriteLengthC1 syn.DistanceC1 ...
                    syn.BranchPointC1 syn.NeuriteLengthC2 ...
                    syn.DistanceC2 syn.BranchPointC2 ...
                    syn.StartXY1 syn.StartXY2];
                tabledata = cat(1,tabledata,row1);
            end
        end
    end
  
    
   end
end

%Utility function
function [p1,p2] = findLinePoints(region)
    rx = region(:,2);
    ry = region(:,1);
    [rxmin,ri1] = min(rx);
    [rxmax,ri2] = max(rx);
    [rymin,ri3] = min(ry);
    [rymax,ri4] = max(ry);
    if (rxmax - rxmin > 0)
        p1 = [rxmin, ry(ri1)];
        p2 = [rxmax, ry(ri2)];
    else
        p1 = [rx(ri3), rymin ];
        p2 = [rx(ri4), rymax ];
    end
end

%Find length of region - assumes straight line - TODO - curve line
function c = findLength(region, resolution)
    sortedregion = sortrows(region,1);
    top = sortedregion(1,:);
    bottom = sortedregion(end,:);
    x1 = top(1);
    y1 = top(2);
    x2 = bottom(1);
    y2 = bottom(2);
    c = sqrt((x2 - x1)^2 + (y2-y1)^2);
    c = c * resolution;
end

%if x,y is at either end of segment within nrange of pixels as [1 1]
function a = isEndpoint(x,y, segment, nrange)
    [m,mrc] = min(segment); %m is min values in each col, mrc is [row,col] of min val
    [n,nrc] = max(segment);
    min1 = segment(mrc(1),:);
    max1 = segment(nrc(1),:);
    a= 0;
    b = [y,x];
      
    if (length(b(abs(b-min1)<= nrange))==2 || length(b(abs(b-max1)<= nrange))==2)
        a = 1;
    end
        
end

function [a12,a13,a23] = findTriangleAngles(p1,p2,p3)
    %extract lengths of triangle sides
    s12 = norm(p2 - p1);
    s13 = norm(p3 - p1);
    s23 = norm(p3 - p2);
    % theta = arccos((a^2 + b^2 - c^2)/2*a*b)
    a12 = acosd((s13^2 + s23^2 - s12^2)/(2*s13*s23))
    a13 = acosd((s12^2 + s23^2 - s13^2)/(2*s12*s23))
    a23 = acosd((s12^2 + s13^2 - s23^2)/(2*s12*s13))
    
end

% Returns row index matching xy coords
function fR=findCSVIndex(xC,yC,StartX,StartY,tol)
    
    d = 0;
    fR = find((StartX >= xC-tol) & (StartX <=xC+tol))
    %if only one val
    if (isempty(fR))
        fR = find((StartY >= yC-tol) & (StartY <=yC+tol))
    end
    
    if (~isempty(fR) && (length(fR) > 1))
        P = [xC,yC]; 
        T = [StartX(fR) StartY(fR)];
        u=tol+1;
        s=0;
        for i=1:length(T)
            X = cat(1,P,T(i,:))
            d = pdist(X,'euclidean')
            if (d <= tol && u > d)
                    u = d;
                    s = fR(i);
            end
        end
        if (s > 0)
            fR =s
        end
    end
    %If still empty or multiple
     if (isempty(fR) || length(fR) > 1)
         fR = [];
     end
    
end
