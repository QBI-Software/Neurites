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
         obj.maskedR(~obj.Iroi) = 0;  % Set all non-keeper pixels to zero.

         obj.maskedG = obj.gbw; % Simply a copy at first.
         obj.maskedG(~obj.Iroi) = 0;  % Set all non-keeper pixels to zero.
         
      end
      
      function obj = set.Segments(obj,segments)
          if (~isempty(segments))
             obj.Segments = segments;
          else
             error('Empty set')
          end
      end
      function obj = set.Synapses(obj,segments)
          if (~isempty(segments))
             obj.Synapses = segments;
          else
             error('Empty set')
          end
      end
      function obj = findSegments1(obj,minlength,tol)
          Rb = bwboundaries(obj.maskedR);
          Gb = bwboundaries(obj.maskedG);
          n1 = size(Rb, 1);
          n2 = size(Gb, 1);
          match = 1;
          segments = {};
          minsynapselength=4;
          for i = 1:n1
            redSegment = Rb{i};
            if (length(redSegment) > minlength)
                x1 = redSegment(:,2);
                y1 = redSegment(:,1);
                for k = 1:n2
                    greenSegment = Gb{k};
                    if (length(greenSegment) <= minlength)
                        continue
                    end
                    x2 = greenSegment(:,2);
                    y2 = greenSegment(:,1);
                    nx1 = size(x1,1);
                    nx2 = size(x2,1);

                    % calculate all pairwise distances
                    %find gives indices which match values
                    K = find(x2 > (min(x1) - tol) & x2 < (max(x1) + tol));
                    O = find(y2 > (min(y1) - tol) & y2 < (max(y1) + tol));
                    %find matching indices in x2,y2
                    if (~isempty(K) && ~isempty(O))
                        changes = zeros(nx1,nx2);
                        for m = 1:nx2
                            %changes(:,m) = sqrt((x1-x2(m)).^2 + (y1-y2(m)).^2);
                            changes(:,m) = (abs(x1-x2(m)) + abs(y1-y2(m)))/2;
                        end
                        [j,n] = find(changes <= tol);
                        if (~isempty(j) && length(j) > minsynapselength)

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
                                    
                                    s  = regionprops(L1, 'centroid');
                                    if (~inpolygon(s.Centroid(1),s.Centroid(2),redSegment(:,2),redSegment(:,1)))
                                        continue
                                    end
                                    if toggle
                                        segments{match}{1} = redSegment;
                                        segments{match}{2} = greenSegment;
                                        segments{match}{3} = [col,row];
                                        toggle = ~toggle;
                                    else
                                        segments{match}{4} = [col,row];
                                        toggle = ~toggle;
                                        match = match+1;
                                    end
                                end
                            else

                                segments{match}{1} = redSegment;
                                segments{match}{2} = greenSegment;
                                segments{match}{3} = [y1(j),x1(j)];
                                segments{match}{4} = [y2(n),x2(n)];
                                match = match + 1;
                            end
                            
                        end
                        %one match per redsegment
                        continue
                    end
                end
            end
          end
          sprintf('Matched=%d of %d', match, n2)
          obj.Segments = segments;
      end
    function obj = findSegments(obj,minlength,tol)
          Rb = bwmorph('skel',obj.rbw);
          Gb = bwmorph(obj.gbw);
          n1 = size(Rb, 1);
          n2 = size(Gb, 1);
          match = 1;
          segments = {};
          minsynapselength=4;
          for i = 1:n1
            redSegment = Rb{i};
            if (length(redSegment) > minlength)
                x1 = redSegment(:,2);
                y1 = redSegment(:,1);
                for k = 1:n2
                    greenSegment = Gb{k};
                    if (length(greenSegment) <= minlength)
                        continue
                    end
                    x2 = greenSegment(:,2);
                    y2 = greenSegment(:,1);
                    nx1 = size(x1,1);
                    nx2 = size(x2,1);

                    % calculate all pairwise distances
                    %find gives indices which match values
                    K = find(x2 > (min(x1) - tol) & x2 < (max(x1) + tol));
                    O = find(y2 > (min(y1) - tol) & y2 < (max(y1) + tol));
                    %find matching indices in x2,y2
                    if (~isempty(K) && ~isempty(O))
                        changes = zeros(nx1,nx2);
                        for m = 1:nx2
                            %changes(:,m) = sqrt((x1-x2(m)).^2 + (y1-y2(m)).^2);
                            changes(:,m) = (abs(x1-x2(m)) + abs(y1-y2(m)))/2;
                        end
                        [j,n] = find(changes <= tol);
                        if (~isempty(j) && length(j) > minsynapselength)

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
                                    
                                    s  = regionprops(L1, 'centroid');
                                    if (~inpolygon(s.Centroid(1),s.Centroid(2),redSegment(:,2),redSegment(:,1)))
                                        continue
                                    end
                                    if toggle
                                        segments{match}{1} = redSegment;
                                        segments{match}{2} = greenSegment;
                                        segments{match}{3} = [col,row];
                                        toggle = ~toggle;
                                    else
                                        segments{match}{4} = [col,row];
                                        toggle = ~toggle;
                                        match = match+1;
                                    end
                                end
                            else

                                segments{match}{1} = redSegment;
                                segments{match}{2} = greenSegment;
                                segments{match}{3} = [y1(j),x1(j)];
                                segments{match}{4} = [y2(n),x2(n)];
                                match = match + 1;
                            end
                            
                        end
                        %one match per redsegment
                        continue
                    end
                end
            end
          end
          sprintf('Matched=%d of %d', match, n2)
          obj.Segments = segments;
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
            if (isequaln(a,b) || isequaln(a,c) || length(a(isnan(a))) > 0)
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
    
    function obj = measureSynapses(obj,showtypes,cell1file,cell2file)
        %TODO
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
