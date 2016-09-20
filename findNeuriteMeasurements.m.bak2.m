classdef NeuritesStimulusRegion
    % ROI region of light stimulation
    %   Segment of annulus
    
    properties
        id
        midline         %midline of stim region (degrees)
        slength         %length of stim region arc (degrees)
        sarea           %area of stim region (px2)
        narea           %area of neurites within region (px2)
        bwroi           %binary roi image
        neurites        %tree, branch, length, area - for dendrogram
        color           %plot colour
        scale
        shiftx
        shifty
        type
    end
    
    methods
        function obj = NeuritesStimulusRegion(id, midline, sarea, slength, neuritesarea, bwroi, color, type)
            obj.id = id;
            obj.midline = midline;
            obj.slength = slength;
            obj.sarea = sarea;
            obj.narea = neuritesarea;
            obj.bwroi = bwroi;
            obj.color = color;  
            obj.type = type;
        end
        function obj = setScale(obj, scale,shiftx,shifty)
            obj.scale = scale;
            obj.shiftx = shiftx;
            obj.shifty = shifty;
        end
        
        %Pass in CSVfile for lookup  
        function obj = analyseROI(obj, scale,shiftx,shifty,CSVfile,neuron)
            fit = 10;
            idx = 1;
            obj = obj.setScale(scale,shiftx,shifty);
            cc = bwconncomp(obj.bwroi)
            [B,L,N,A] = bwboundaries(obj.bwroi,4,'noholes');
            for k = 1:length(B)
                boundary = B{k};
                plot(boundary(:,2), boundary(:,1), obj.color, 'LineWidth', 1);
                
                [xC, yC]= obj.img2Coords(s.Extrema(k,1),s.Extrema(k,2));
            end
            
            figure
            plot(stats(1).PixelList, 'xr')
            frcache = containers.Map();
            for k = 1:length(stats(1).PixelList)
                [xC, yC]= obj.img2Coords(stats(1).PixelList(k,1),stats(1).PixelList(k,2));
            
%             s = regionprops(obj.bwroi,'Extrema','BoundingBox');
%             xb = [s.BoundingBox(1),s.BoundingBox(1)+s.BoundingBox(3),s.BoundingBox(1)+s.BoundingBox(3),s.BoundingBox(1),s.BoundingBox(1)]
%             yb = [s.BoundingBox(2),s.BoundingBox(2),s.BoundingBox(2)+s.BoundingBox(4),s.BoundingBox(2)+s.BoundingBox(4),s.BoundingBox(2)];
%             hold on
%             plot(s.Extrema(:,1),s.Extrema(:,2),'xr','LineWidth',2)
%             plot(xb,yb,'-g', 'LineWidth', 3)
%             %Convert boundingbox
%             xb = (xb-shiftx)/scale;
%             yb = -(yb+shifty)/scale;
%             frcache = containers.Map();

%             for k = 1:length(s.Extrema)
                [xC, yC]= obj.img2Coords(s.Extrema(k,1),s.Extrema(k,2));
                %May need to adjust fit param to find CSV data (depends on how well overlaid
                %Load to table? 
                
                iterate = 5; 
                for j = 1:iterate
                    fR=findCSVIndex(xC,yC,CSVfile.StartX,CSVfile.StartY,fit);
                    if (isempty(fR))
                        fit = fit + 2;
                    else
                        a = num2str(fR);
                        if (~frcache.isKey(a))%empty(~find(ismember(frcache,fR))))
                            %Identify node               
                            t = CSVfile.Tree(fR)
                            b = CSVfile.Order(fR)
                            node = neuron{t,1}{1,b}
                            for i=1:length(node)
                                n = node(i,1)
                                p = pdist2([xC yC],n.points)
                                [r,c] = find(p < 10) %C IS IDX OF N.points matching
                                if (length(c)>0)
                                    if(length(c)>1)
                                        [ei,e] = min(p(c))
                                        c=e
                                    end
                                    matchp = c %match coords
                                    startxy = n.points(c,:)
                                    endxy = n.points(c,:);
                                    %find startxy
                                    c=c-1;
                                    while (c > 0 & inpolygon(n.points(c,1),n.points(c,2),xb',yb'))
                                        startxy = n.points(c,:)
                                        c =c-1;
                                    end
                                    
                                    %find endxy
                                    c = matchp + 1
                                    while (c <= length(n.points) & inpolygon(n.points(c,1),n.points(c,2),xb',yb'))
                                        endxy = n.points(c,:)
                                        c =c+1;
                                    end
                                    %Test view
                                    figure
                                    hold on
                                    plot(xb',yb','-g', 'LineWidth', 3)
                                    plot(startxy(1),startxy(2),'xr')
                                    plot(endxy(1),endxy(2),'or')
                                    plot(n.points(:,1), n.points(:,2),'-c');
                                    %Calculations for segment
                                    [len,vol,area,endpoints,branchpoints] = findNeuriteMeasurements(obj,startxy(1),startxy(2),endxy,CSVfile,fR);
                                    [dsoma,vsoma,ssoma,somapoints] = findSomaMeasurePoints(obj,xC,yC,CSVfile,fR);
                                    nlength = pdist2(startxy,endxy) %OR cumulative loop
                                    marker= struct('fr',fR,'x',xC,'y',yC,'nodeid',n.id,'tree',t,'branch',b,'blength',n.branchlength,'nlength',len)
                                    frcache(a) = marker; % = cat(1,frcache,fR);
                                    obj.neurites(idx,1).crow = fR;
                                    obj.neurites(idx,1).tree = t;
                                    obj.neurites(idx,1).branch = b;
                                    obj.neurites(idx,1).blength = n.branchlength;
                                    obj.neurites(idx,1).nlength = len;
                                    %obj.neurites(idx,1).nvol = vol;
                                    %obj.neurites(idx,1).nsa = area;
                                    obj.neurites(idx,1).npoints = endpoints;        %points to end of neurite for overlay
                                    obj.neurites(idx,1).xy = [xC,yC];              %searched coordinates
                                    obj.neurites(idx,1).somad = dsoma;
                                    %obj.neurites(idx,1).somav = vsoma;
                                    %obj.neurites(idx,1).somas = ssoma;
                                    obj.neurites(idx,1).somapoints = somapoints;    %points back to soma for overlay
                                    idx = idx + 1;

                                end
                            end
                        end
                    break %break from fR iteration    
                    end
                end
            end
            
        end
        %Requires setting  scale,shiftx,shifty,
        function [x,y] = img2Coords(obj,imgx, imgy)
            x = (imgx - obj.shiftx)/obj.scale;
            y = (-imgy - obj.shifty)/obj.scale;

          end
          function [imgx,imgy] = coords2Img(obj,x, y)
            imgx = (x * obj.scale) + obj.shiftx;
            imgy = -((y * obj.scale) + obj.shifty);

          end
    end
    
end


