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
             error('Error: Cannot set SEGMENTS - Empty set')
          end
      end
      function obj = set.Synapses(obj,segments)
          if (~isempty(segments))
             obj.Synapses = segments;
          else
             error('Error: Cannot set SYNAPSES - Empty set')
          end
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
          
          %test
          %figure
          imshow(obj.maskedR)
          hold on
          %Masked region coords
          [mar,mbr] = find(obj.maskedR==1);
          [mag,mbg] = find(obj.maskedG==1);
          for i = 1:nR
            redSegment = Rb{i};
            %Check is in ROI
            if (isempty(find(ismember(redSegment,[mar,mbr])==1, 1)))
                %plot(redSegment(:,2),redSegment(:,1),'b','LineWidth',2)
                continue
            else
                plot(redSegment(:,2),redSegment(:,1),'r','LineWidth',2)
            end
            if (length(redSegment) > minlength)
                x1 = redSegment(:,2);
                y1 = redSegment(:,1);
                for k = 1:nG
                    greenSegment = Gb{k};
                    %Check valid ROI
                    if (isempty(find(ismember(greenSegment,[mag,mbg])==1,1)))
                       % plot(greenSegment(:,2),greenSegment(:,1),'y','LineWidth',2)
                        continue
                    end
                    x2 = greenSegment(:,2);
                    y2 = greenSegment(:,1);
                    nx1 = size(x1,1);
                    nx2 = size(x2,1);

                    
                    
                    % Filter first - gives indices which match values
                    Gxi = find(x2 >= (min(x1) - tol) & x2 <= (max(x1) + tol));
                    Gyi = find(y2 >= (min(y1) - tol) & y2 <= (max(y1) + tol));
                    [~,Gi] = intersect(Gxi,Gyi);
                    %find matching indices in x2,y2
                    if (~isempty(Gi))
                        tic
                        changes = zeros(nx1,nx2);
                        for m = 1:nx2
                            changes(:,m) = (abs(x1-x2(m)) + abs(y1-y2(m)))/2;
                        end
                        [rr,gc] = find(changes <= tol); % idx = red row , green col
                        if (~isempty(rr))

                            %check that range is single line
                            blankimg = zeros(size(obj.rbw));
                            blankimg(y1(rr),x1(rr)) = 1;
                            blankimg(y2(gc),x2(gc)) = 1;
                            [L, num] = bwlabel(blankimg, 8);
                            for w=1:num
                                % DETECT R and G SYNAPTIC REGION
                                [ry,cx] = find(L==w);
                                %L1 = L == w;
                                %[c,r] = find(L1);
                                rc = [ry,cx];
                                arr=find(ismember(rc,[y1(rr),x1(rr)],'rows'));
                                agg=find(ismember(rc,[y2(gc),x2(gc)],'rows'));
                                if (~isempty(arr) && ~isempty(agg))
                                    segments{match}{1} = redSegment;
                                    segments{match}{2} = greenSegment;
                                    segments{match}{3} = rc(arr,:);
                                    segments{match}{4} = rc(agg,:);
                                    segments{match}{5} = findLength(rc(arr,:), resolution);
                                    segments{match}{6} = findLength(rc(agg,:), resolution);
                                    match = match + 1;
                                end
                            end
                            
                        end
                        toc
                        %one match per redsegment
                        %continue
                    end
                end %greenseg loop
            end
          end %redseg loop
          
          hold off
          if (~isempty(segments))
              obj.Segments = segments;
              status = sprintf('Found %d Matches', match)
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
           % Show on image
            plot(redSegment(:,2), redSegment(:,1), 'r', 'LineWidth', 2);
            plot(greenSegment(:,2), greenSegment(:,1),'g','LineWidth',2);
            plot(redhighlight(:,2), redhighlight(:,1), 'y', 'LineWidth', 3);
            plot(greenhighlight(:,2), greenhighlight(:,1),'y','LineWidth',3);

            %find min max points on each line
            [p1,p2] = findLinePoints(redhighlight);
            [p3,p4] = findLinePoints(greenhighlight);

            
            %determine median synapse location
            rmarkerx = median(redhighlight(:,2));
            rmarkery = median(redhighlight(:,1));
            gmarkerx = median(greenhighlight(:,2));
            gmarkery = median(greenhighlight(:,1));
%             nrange = [1 1];
%             %endpoint of one cell on another
%             if (isEndpoint(rmarkerx, rmarkery, redSegment,nrange) || ...
%                 isEndpoint(gmarkerx,gmarkery,greenSegment,nrange))
%                 type = 2;
%             end
            mycolors = ['m','b','c'];
            if (showplots && ismember(type,showtypes))
                plot(p1(1), p1(2), p2(1), p2(2), p3(1), p3(2), p4(1), p4(2),...
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
%         I = obj.rbw;
        T1 = readtable(cell1file);
        T2 = readtable(cell2file);
        %plot overlay
%         X = (T1.StartX * scale) + shiftx;
%         Y = (T1.StartY * scale) + shifty;
%         figure
%         imshow(I);
%         hold on;
%         plot(X,-Y,'color','c','LineStyle','-','LineWidth', 2);
%         hold off;
        %tol = 10;
        saved ={};
        ctr=1;
        for i = 1: length(obj.Synapses)
            syn = obj.Synapses{i};
            
            if (ismember(syn.SynapseType,showtypes))
                %cell1
                syn.shiftx = shiftx;
                syn.shifty = shifty;
                syn.scale = scale;
                
                %median coords - Cell1
                [xC, yC]= syn.img2Coords(syn.MedianC1(1),syn.MedianC1(2));
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
                    syn.NeuriteLengthC1 = T1.Length__m_(fR);
                    syn.DistanceC1= findDistanceToSoma(xC,yC,syn.StartXY1,T1.LengthToBeginning__m_(fR));
                    syn.BranchPointC1 = T1.Order(fR);
                    syn.BranchTypeC1 = T1.PointType(fR);
                    [syn.NeuriteEndC1,syn.EndpointsC1] = findLength2NeuriteEnd(xC,yC,T1,fR);
                end
                
                %median coords - Cell2
                [xC, yC]= syn.img2Coords(syn.MedianC2(1),syn.MedianC2(2));
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
                    syn.NeuriteLengthC2 = T2.Length__m_(fR);
                    syn.DistanceC2= findDistanceToSoma(xC,yC,syn.StartXY2,T2.LengthToBeginning__m_(fR));
                    syn.BranchPointC2 = T2.Order(fR);
                    syn.BranchTypeC2 = T2.PointType(fR);
                    [syn.NeuriteEndC2,syn.EndpointsC2] = findLength2NeuriteEnd(xC,yC,T2,fR);
                end
                %check if duplicate
                syn
                if (~obj.duplicateSynapse(syn, i))
                    %save it
                    %obj.Synapses{i} = syn;
                    saved{ctr} = syn;
                    ctr = ctr + 1;
                end
                
            end
        end
       
    %length(obj.Synapses)
    obj.Synapses = saved;
    length(obj.Synapses)

    end
    function [colnames,tabledata] = generateTable(obj,types, cell1label, cell2label)
         colnames = {'Type' 'Cell1_X' 'Cell1_Y' 'Cell2_X' 'Cell2_Y',...
             'Cell1_length' 'Cell1_distance' 'Cell1_order' 'Cell1_end' ...
             'Cell2_length' 'Cell2_distance' 'Cell2_order' 'Cell2_end' ...
             'Cell1_StartX' 'Cell1_StartY' 'Cell2_StartX' 'Cell2_StartY' };
         colnames = strrep(colnames, 'Cell1', cell1label);
         colnames = strrep(colnames, 'Cell2', cell2label);
%         colnames = {'Type' 'DS_X' 'DS_Y' 'SBAC_X' 'SBAC_Y' ...
%             'DS_length' 'DS_distance' 'DS_order' 'DS_end'  ...
%             'SBAC_length' 'SBAC_distance' 'SBAC_order' 'SBAC_end' ...
%             'DS_StartX' 'DS_StartY' 'SBAC_StartX' 'SBAC_StartY' };
        tabledata =[];
        for i=1:length(obj.Synapses)
            syn = obj.Synapses{i};
            if (ismember(syn.SynapseType,types))
                
                %show data in table
                row1 = [syn.SynapseType ...
                    syn.MedianC1 syn.MedianC2 ...
                    syn.NeuriteLengthC1 syn.DistanceC1 ...
                    syn.BranchPointC1 syn.NeuriteEndC1 syn.NeuriteLengthC2 ...
                    syn.DistanceC2 syn.BranchPointC2 syn.NeuriteEndC2 ...
                    syn.StartXY1 syn.StartXY2];
                tabledata = cat(1,tabledata,row1);
            end
        end
    end
    %Check if already loaded - num: is iteration of updated Syn objects
    function isduplicate = duplicateSynapse(obj,syn, num)
        isduplicate = 0;
        for i=1:num
            synA = obj.Synapses{i};
            if (synA.isequal(syn))
                isduplicate = 1;
                continue;
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
    if(isempty(region))
        c = 0;
    else
        sortedregion = sortrows(region,1);
        top = sortedregion(1,:);
        bottom = sortedregion(end,:);
        x1 = top(1);
        y1 = top(2);
        x2 = bottom(1);
        y2 = bottom(2);
        V = [x1 y1;x2 y2];
        c = pdist(V,'euclidean') * resolution;
    end
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
    fR = find((StartX >= xC-tol) & (StartX <=xC+tol));
    %if only one val
    if (isempty(fR))
        fR = find((StartY >= yC-tol) & (StartY <=yC+tol));
    end
    
    if (~isempty(fR) && (length(fR) > 1))
        P = [xC,yC]; 
        T = [StartX(fR) StartY(fR)];
        u=tol+1;
        s=0;
        for i=1:length(T)
            X = cat(1,P,T(i,:));
            d = pdist(X,'euclidean')
            if (d <= tol && u > d)
                    u = d;
                    s = fR(i);
            end
        end
        if (s > 0)
            fR =s;
        end
    end
    %If still empty or multiple
     if (isempty(fR) || length(fR) > 1)
         fR = [];
     end
    
end

function d = distanceBetween(A, B, type)
    d= 0;
    if (isempty(type))
        type = 'euclidean';
    end
    X = cat(1,A,B);
    d = pdist(X,type);
end

function d = findDistanceToSoma(x,y,startxy,lengthxy)
    d = distanceBetween([x,y],startxy,'euclidean');
    d = lengthxy - d;
end

function [d,endpoints] = findLength2NeuriteEnd(x,y,T,fR)
    maxfr = height(T);
    endxy = [T.EndX(fR), T.EndY(fR)];
    synxy = [x,y];
    %endpoints = {}
    xpoints = [T.StartX(fR)];
    ypoints = [T.StartY(fR)];
    d = distanceBetween(synxy,endxy,'euclidean');
    while(strcmp(T.PointType(fR),'EP')== 0 && fR < maxfr)
        d = d + T.Length__m_(fR);
        %Save points for display
        xpoints(end+1) = T.StartX(fR);
        ypoints(end+1) = T.StartY(fR);
        fR = fR + 1;
    end    
    endpoints =[xpoints,ypoints];
end

      
