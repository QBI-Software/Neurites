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
        neuron          %digital neuron model nodes
    end
    
    methods
        function obj = NeuritesStimulusRegion(id, midline, sarea, slength, neuritesarea, bwroi, color, type,neuron)
            obj.id = id;
            obj.midline = midline;
            obj.slength = slength;
            obj.sarea = sarea;
            obj.narea = neuritesarea;
            obj.bwroi = bwroi;
            obj.color = color;  
            obj.type = type;
            obj.neuron = neuron;
        end
        function obj = setScale(obj, scale,shiftx,shifty)
            obj.scale = scale;
            obj.shiftx = shiftx;
            obj.shifty = shifty;
        end
        
        %Pass in CSVfile for lookup  
        function obj = analyseROI(obj, CSVfile,FR,dcache)
            %fit = 10;
            idx = 1;
            %obj = obj.setScale(scale,shiftx,shifty);
            %frcache = containers.Map();
            
            k=1;
            while k<=length(FR)
                fR = FR(k);
                obj.neurites(idx,1).crow = fR;
                xC = CSVfile.StartX(fR);
                yC = CSVfile.StartY(fR);
                prelength = 0;
                if (fR <= length(dcache) && isfield(dcache(fR),'begin') && ~isempty(dcache(fR).begin))
                    [beginx,beginy] = obj.img2Coords(dcache(fR).begin(1),dcache(fR).begin(2));
                    prelength = pdist2([beginx,beginy],[xC,yC]);
                end
                dsoma = CSVfile.LengthToBeginning__m_(fR);
                %[dsoma,vsoma,ssoma,somapoints] = findSomaMeasurePoints(obj,xC,yC,CSVfile,fR);
                dsoma = dsoma - prelength
                tree = CSVfile.Tree(fR)
                branch = CSVfile.Order(fR)
                branchlength = obj.findBranchlength(tree,branch,xC,yC)
                len = 0;
                %find end xy 
                while(fR<max(FR) && tree == CSVfile.Tree(fR) && branch== CSVfile.Order(fR))
                    len = len + CSVfile.Length__m_(fR);
                    k = k+1;
                    fR = FR(k);
                end
                k = k+1;
                postlength=0;
                if (fR <= length(dcache) && isfield(dcache(fR),'end') && ~isempty(dcache(fR).end))
                    [endx,endy] = obj.img2Coords(dcache(fR).end(1),dcache(fR).end(2));  
                    postlength = pdist2([endx,endy],[CSVfile.EndX(fR),CSVfile.EndY(fR)]);
                end
                %[len,vol,area,endpoints,branchpoints] = findNeuriteMeasurements(obj,xC,yC,endxy,CSVfile,fR);
                len = len + prelength + postlength
                %obj.neurites(idx,1).crow = fR;
                obj.neurites(idx,1).tree = tree;
                obj.neurites(idx,1).branch = branch;
                obj.neurites(idx,1).blength = branchlength;
                obj.neurites(idx,1).nlength = len;
                %obj.neurites(idx,1).nvol = vol;
                %obj.neurites(idx,1).nsa = area;
                %obj.neurites(idx,1).npoints = endpoints;        %points to end of neurite for overlay
                obj.neurites(idx,1).xy = [xC,yC];              %searched coordinates
                obj.neurites(idx,1).somad = dsoma;
                %obj.neurites(idx,1).somav = vsoma;
                %obj.neurites(idx,1).somas = ssoma;
                %obj.neurites(idx,1).somapoints = somapoints;    %points back to soma for overlay
                idx = idx + 1;
                
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
          
          function branchlength = findBranchlength(obj,tree,branch,xC,yC)
            branchlength = 0;
            if ~isempty(obj.neuron)
              %Identify node  
                node = obj.neuron{tree,1}{1,branch};
                for i=1:length(node)
                    n = node(i,1)
                    p = pdist2([xC yC],n.points);
                    [r,c] = find(p < 10); %C IS IDX OF N.points matching
                    if(~isempty(c))
                        branchlength = n.branchlength;
                        break
                    end
                end
            else
                disp('Cannot determine branchlength - neuron not loaded')
            end
          end
    end
    
end


