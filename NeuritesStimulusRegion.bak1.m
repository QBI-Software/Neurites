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
        function obj = analyseROI(obj, scale,shiftx,shifty,CSVfile,neuron,FR)
            fit = 10;
            idx = 1;
            obj = obj.setScale(scale,shiftx,shifty);
            frcache = containers.Map();
            
            
            
            
%             cc = bwconncomp(obj.bwroi) %number connected objs
%             [B,L,N,A] = bwboundaries(obj.bwroi,4,'noholes');
%              [BY,BX] = find(obj.bwroi ==1);
%              plot(BX(:,1),BY(:,1),'xr') %check
%              xb = (BX(:,1)-shiftx)/scale;
%              yb = -(BY(:,1)+shifty)/scale;
            
            tree = 0;
            branch = 0;
            pts =[];
            iterate = 5; 
           
            for k=1:length(FR)
                fR = FR(k);
                xC = CSVfile.StartX(fR);
                yC = CSVfile.StartY(fR);
                tree = CSVfile.Tree(fR);
                branch = CSVfile.Order(fR);
                %find end xy 
                while(tree == CSVfile.Tree(fR) && branch== CSVfile.Order(fR))
                    k = k+1;
                    fR = FR(k);
                end
                k = k-1;
                fR = FR(k);
                endxy = [CSVfile.EndX(fR),CSVfile.EndY(fR)]
                [len,vol,area,endpoints,branchpoints] = findNeuriteMeasurements(obj,xC,yC,endxy,CSVfile,fR);
                [dsoma,vsoma,ssoma,somapoints] = findSomaMeasurePoints(obj,xC,yC,CSVfile,fR);
                obj.neurites(idx,1).crow = fR;
                obj.neurites(idx,1).tree = tree;
                obj.neurites(idx,1).branch = branch;
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
                %Identify node  
                node = neuron{tree,1}{1,branch};
                for i=1:length(node)
                    n = node(i,1)
                    p = pdist2([xC yC],n.points)
                    [r,c] = find(p < 10) %C IS IDX OF N.points matching
                    if(length(c)>0)
                        break
                    end
                end
                if ~frcache.isKey(n.id)
                    frcache(n.id) = n;
                    
                    len = len + CSVfile.Length__m_(fR)
                    somad = CSVfile.LengthToBeginning__m_(fR)
                    [len,vol,area,endpoints,branchpoints] = findNeuriteMeasurements(obj,startxy(1),startxy(2),endxy,CSVfile,fR);
                    [dsoma,vsoma,ssoma,somapoints] = findSomaMeasurePoints(obj,xC,yC,CSVfile,fR);

                                    %marker= struct('fr',fR,'x',xC,'y',yC,'nodeid',n.id,'tree',t,'branch',b,'blength',n.branchlength,'nlength',len)
                                    %frcache(a) = marker; % = cat(1,frcache,fR);
                                    obj.neurites(idx,1).crow = fR;
                                    obj.neurites(idx,1).tree = tree;
                                    obj.neurites(idx,1).branch = branch;
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

                                    %Find corresponding point in boundary points
                                    p = pdist2(endxy, [xb,yb]);
                                    [val,idx] = find(p < 10)
                                    if(length(idx)>1)
                                        [ei,k] = min(p(idx)) % sets next boundary idx
                                    end
                                end
                            end
                        end
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


