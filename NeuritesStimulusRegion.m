classdef NeuritesStimulusRegion
    % ROI region of light stimulation
    %   Segment of annulus
    
    properties
        midline         %midline of stim region (degrees)
        slength         %length of stim region arc (degrees)
        sarea           %area of stim region (px2)
        narea           %area of neurites within region (px2)
        boundaries      %boundaries from bwboundaries - may be individual or joined neurites
        neurites        %tree, branch, length - for dendrogram
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
            obj = obj.setScale(scale,shiftx,shifty);
            for k = 1:length(obj.boundaries)
                boundary = obj.boundaries{k};
                %find min max points on each line
                [minp,maxp] = findLinePoints(boundary);
                %Check
                hold on
                plot(minp(1),minp(2),'xr')
                plot(maxp(1),maxp(2),'xr')
                %Get data for min
                [xC, yC]= obj.img2Coords(minp(1),minp(2));
                fR=findCSVIndex(xC,yC,CSVfile.StartX,CSVfile.StartY,fit);
                
                obj.neurites(k,1).tree = CSVfile.Tree(fR);
                obj.neurites(k,1).branch = CSVfile.Order(fR);
                obj.neurites(k,1).nlength = CSVfile.Length__m_(fR);
                obj.neurites(k,1).somax = findDistanceToSoma(xC,yC,[CSVfile.StartX(fR), CSVfile.StartY(fR)], CSVfile.LengthToBeginning__m_(fR));
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


