classdef NeuritesStimulusRegion
    % ROI region of light stimulation
    %   Segment of annulus
    
    properties
        midline         %midline of stim region (degrees)
        slength         %length of stim region arc (degrees)
        sarea           %area of stim region (px2)
        narea           %area of neurites within region (px2)
        boundaries      %boundaries from bwboundaries - may be individual or joined neurites
        neurites        %tree, branch, length, area - for dendrogram
        color           %plot colour
        scale
        shiftx
        shifty
    end
    
    methods
        function obj = NeuritesStimulusRegion(midline, sarea, slength, neuritesarea, boundaries, color)
            obj.midline = midline;
            obj.slength = slength;
            obj.sarea = sarea;
            obj.narea = neuritesarea;
            obj.boundaries = boundaries;
            obj.color = color;       
        end
        function obj = setScale(obj, scale,shiftx,shifty)
            obj.scale = scale;
            obj.shiftx = shiftx;
            obj.shifty = shifty;
        end
        
        %Pass in CSVfile for lookup  
        function obj = analyseBoundaries(obj, scale,shiftx,shifty,CSVfile)
            fit = 10;
            idx = 1;
            obj = obj.setScale(scale,shiftx,shifty);
            for k = 1:length(obj.boundaries)
                boundary = obj.boundaries{k};
                %find min max points on each line
                [minp,maxp] = findLinePoints(boundary);
                %Check
                hold on
                plot(minp(1),minp(2),'xr')
                plot(maxp(1),maxp(2),'xy')
                
                %Get data for min
                [xC, yC]= obj.img2Coords(minp(1),minp(2));
                iterate = 5; %May need to adjust fit param to find CSV data (depends on how well overlaid)
                for j = 1:iterate
                    fR=findCSVIndex(xC,yC,CSVfile.StartX,CSVfile.StartY,fit)
                    if (isempty(fR))
                        fit = fit + 2;
                    else
                        break
                    end
                end
                if (~isempty(fR))
                    obj.neurites(idx,1).crow = fR;
                    obj.neurites(idx,1).tree = CSVfile.Tree(fR);
                    obj.neurites(idx,1).branch = CSVfile.Order(fR);
                    obj.neurites(idx,1).nlength = CSVfile.Length__m_(fR);
                    obj.neurites(idx,1).boundary = boundary; %entire branch pixel area
                    %obj.neurites(idx,1).somax = findDistanceToSoma(xC,yC,[CSVfile.StartX(fR), CSVfile.StartY(fR)], CSVfile.LengthToBeginning__m_(fR));
                    [d,v,s,points] = findSomaMeasurePoints(obj,xC,yC,CSVfile,fR);
                    obj.neurites(idx,1).somad = d;
                    obj.neurites(idx,1).somav = v;
                    obj.neurites(idx,1).somas = s;
                    obj.neurites(idx,1).somapoints = points; %points back to soma for overlay
                end
                %Get data for max - save if different branch
                [xC, yC]= obj.img2Coords(maxp(1),maxp(2));
                for j = 1:iterate
                    fRmax=findCSVIndex(xC,yC,CSVfile.StartX,CSVfile.StartY,fit)
                    if (isempty(fRmax))
                        fit = fit + 2;
                    else
                        break
                    end
                end
                
                if (~isempty(fRmax) && fRmax ~= fR)
                    disp('Different row detected for max');
                    %idx = idx+1;
                    fR = fRmax;
                    %somax = findDistanceToSoma(xC,yC,[CSVfile.StartX(fR), CSVfile.StartY(fR)], CSVfile.LengthToBeginning__m_(fR));
                    [d,v,s,points] = findSomaMeasurePoints(obj,xC,yC,CSVfile,fR);
                    %if different tree then add as separate row
                    if (CSVfile.Tree(fR) ~= obj.neurites(idx,1).tree)
                        idx = idx+1;
                        obj.neurites(idx,1).tree = CSVfile.Tree(fR);
                        obj.neurites(idx,1).branch = CSVfile.Order(fR);
                        obj.neurites(idx,1).nlength = CSVfile.Length__m_(fR);
                        obj.neurites(idx,1).narea = length(boundary)/scale; 
                        %obj.neurites(idx,1).somax = somax; 
                        obj.neurites(idx,1).somad = d;
                        obj.neurites(idx,1).somav = v;
                        obj.neurites(idx,1).somas = s;
                        obj.neurites(idx,1).somapoints = points;
                    elseif (obj.neurites(idx,1).somad < d) %replace with max point measurements
                        obj.neurites(idx,1).nlength = d - obj.neurites(idx,1).somad;
                        %obj.neurites(idx,1).somax = somax;
                        obj.neurites(idx,1).somad = d;
                        obj.neurites(idx,1).somav = v;
                        obj.neurites(idx,1).somas = s;
                        obj.neurites(idx,1).somapoints = points;
                        obj.neurites(idx,1).crow = [obj.neurites(idx,1).crow, fRmax];%append new rownumber
                        obj.neurites(idx,1).branch = [obj.neurites(idx,1).branch, CSVfile.Order(fR)]; %append new branch
                    elseif (obj.neurites(idx,1).somad > d) %ie reversed min/max
                        obj.neurites(idx,1).nlength = obj.neurites(idx,1).somad - d;
                        obj.neurites(idx,1).crow = [obj.neurites(idx,1).crow, fRmax];%append new rownumber
                        obj.neurites(idx,1).branch = [obj.neurites(idx,1).branch, CSVfile.Order(fR)]; %append new branch
                    end
                end
                idx = idx+1;
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


