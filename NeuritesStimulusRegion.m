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
        function obj = analyseROI(obj, scale,shiftx,shifty,CSVfile)
            fit = 10;
            idx = 1;
            obj = obj.setScale(scale,shiftx,shifty);
            s = regionprops(obj.bwroi,'Extrema');
            hold on
            plot(s.Extrema(:,1),s.Extrema(:,2),'xr','LineWidth',2)
            frcache = containers.Map();
            T = table;%(fR,Tree,Branch,Start,End,Length)
            FR = [];
            TREE = [];
            BRANCH = [];
            STARTX = [];
            STARTY = [];
            for k = 1:length(s.Extrema)
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
                            marker= struct('fr',fR,'x',xC,'y',yC);
                            frcache(a) = marker; % = cat(1,frcache,fR);
                            FR = cat(1,FR,fR);
                            TREE = cat(1,TREE,CSVfile.Tree(fR));
                            BRANCH = cat(1, BRANCH, CSVfile.Order(fR));
                            STARTX = cat(1, STARTX, CSVfile.StartX(fR));
                            STARTY = cat(1, STARTY, CSVfile.StartY(fR));
                        end
                        break
                    end
                end
            end
            T.FR = FR;
            T.TREE = TREE;
            T.BRANCH = BRANCH;
            T.STARTX = STARTX;
            T.STARTY = STARTY;
            T
            statarry = grpstats(T,{'TREE', 'BRANCH'})
            g = statarry(statarry.GroupCount >1, :)
            for i= 1:height(g)
                T1 = T((T.TREE == g.TREE(i) & T.BRANCH ==g.BRANCH(i)),:)
                r1 = frcache(num2str(T1.FR(1)))
                r2 = frcache(num2str(T1.FR(2)))
                xC = r1.x
                yC = r1.y
                xC2 = r2.x
                yC2 = r2.y
                fR = r1.fr
                [l,o,a,endpoints,branchpoints] = findNeuriteMeasurements(obj,xC,yC,[xC2 yC2],CSVfile,fR);
                [d,v,s,somapoints] = findSomaMeasurePoints(obj,xC,yC,CSVfile,fR);
                obj.neurites(idx,1).crow = fR;
                obj.neurites(idx,1).tree = CSVfile.Tree(fR);
                obj.neurites(idx,1).branch = CSVfile.Order(fR);
                 obj.neurites(idx,1).nlength = l;
                 obj.neurites(idx,1).nvol = o;
                 obj.neurites(idx,1).nsa = a;
%                 obj.neurites(idx,1).npoints = 0;%endpoints;        %points to end of neurite for overlay
                obj.neurites(idx,1).xy = [xC,yC];              %searched coordinates
                obj.neurites(idx,1).somad = d;
                obj.neurites(idx,1).somav = v;
                obj.neurites(idx,1).somas = s;
                obj.neurites(idx,1).somapoints = somapoints;    %points back to soma for overlay
                idx = idx + 1;
            end
                    
%             numneurites=length(frcache)
%             fc = sort(frcache.keys)
           % marker1 = frcache(fc(1))
           % marker2 = frcache(fc(2))
             
            %[l,o,a,endpoints,branchpoints] = findNeuriteRegion(obj,marker1,marker2,CSVfile); %neurite segment only
%                     for i=1:length(branchpoints)
%                         fRi = branchpoints(i);
%                         if (fRi ~= fR)
%                             disp('Different row detected for max');
%                             xC = 0;
%                             yC = 0;
%                             [l,o,a,endpoints,branchpoints] = findNeuriteMeasurements(obj,xC,yC,[xC2 yC2],CSVfile,fR);
%                         end
%            for i=1:length(fc)
%                key = fc(i);
%                
%                marker = frcache(key{1})
%                
%                xC = marker.x;
%                yC = marker.y;
%                fR = marker.fr;
%               
%                tree = CSVfile.Tree(fR)
%                branch = CSVfile.Order(fR)
%                l = 0; %blank for non-segment
%                o = 0;
%                a = 0;
%                %Detect markers on same branch
%                if (i < length(fc))
%                    key2 = fc(i+1);
%                    marker2 = frcache(key2{1})
%                    if(CSVfile.Tree(marker2.fr) == tree &&  CSVfile.Order(marker2.fr) == branch) 
%                        xC2 = marker2.x;
%                        yC2 = marker2.y;
%                        [l,o,a,endpoints,branchpoints] = findNeuriteMeasurements(obj,xC,yC,[xC2 yC2],CSVfile,fR);
%                    end
%                end
%                 [d,v,s,somapoints] = findSomaMeasurePoints(obj,xC,yC,CSVfile,fR);
%                 obj.neurites(idx,1).crow = fR;
%                 obj.neurites(idx,1).tree = CSVfile.Tree(fR);
%                 obj.neurites(idx,1).branch = CSVfile.Order(fR);
%                  obj.neurites(idx,1).nlength = l;
%                  obj.neurites(idx,1).nvol = o;
%                  obj.neurites(idx,1).nsa = a;
% %                 obj.neurites(idx,1).npoints = 0;%endpoints;        %points to end of neurite for overlay
%                 obj.neurites(idx,1).xy = [xC,yC];              %searched coordinates
%                 obj.neurites(idx,1).somad = d;
%                 obj.neurites(idx,1).somav = v;
%                 obj.neurites(idx,1).somas = s;
%                 obj.neurites(idx,1).somapoints = somapoints;    %points back to soma for overlay
%                 idx = idx + 1;
                   % end
              %  end
                %Get data for max - save if different branch
%                 [xC, yC]= obj.img2Coords(maxp(1),maxp(2));
%                 for j = 1:iterate
%                     fRmax=findCSVIndex(xC,yC,CSVfile.StartX,CSVfile.StartY,fit)
%                     if (isempty(fRmax))
%                         fit = fit + 2;
%                     else
%                         break
%                     end
%                 end
%                 
%                 if (~isempty(fRmax) && fRmax ~= fR)
%                     disp('Different row detected for max');
%                     %idx = idx+1;
%                     fR = fRmax;
%                     %somax = findDistanceToSoma(xC,yC,[CSVfile.StartX(fR), CSVfile.StartY(fR)], CSVfile.LengthToBeginning__m_(fR));
%                     [d,v,s,points] = findSomaMeasurePoints(obj,xC,yC,CSVfile,fR);
%                     %if different tree then add as separate row
%                     if (CSVfile.Tree(fR) ~= obj.neurites(idx,1).tree)
%                         idx = idx+1;
%                         obj.neurites(idx,1).tree = CSVfile.Tree(fR);
%                         obj.neurites(idx,1).branch = CSVfile.Order(fR);
%                         obj.neurites(idx,1).nlength = CSVfile.Length__m_(fR);
%                         obj.neurites(idx,1).narea = length(boundary)/scale; 
%                         %obj.neurites(idx,1).somax = somax; 
%                         obj.neurites(idx,1).somad = d;
%                         obj.neurites(idx,1).somav = v;
%                         obj.neurites(idx,1).somas = s;
%                         obj.neurites(idx,1).somapoints = points;
%                     elseif (obj.neurites(idx,1).somad < d) %replace with max point measurements
%                         obj.neurites(idx,1).nlength = d - obj.neurites(idx,1).somad;
%                         %obj.neurites(idx,1).somax = somax;
%                         obj.neurites(idx,1).somad = d;
%                         obj.neurites(idx,1).somav = v;
%                         obj.neurites(idx,1).somas = s;
%                         obj.neurites(idx,1).somapoints = points;
%                         obj.neurites(idx,1).crow = [obj.neurites(idx,1).crow, fRmax];%append new rownumber
%                         obj.neurites(idx,1).branch = [obj.neurites(idx,1).branch, CSVfile.Order(fR)]; %append new branch
%                     elseif (obj.neurites(idx,1).somad > d) %ie reversed min/max
%                         obj.neurites(idx,1).nlength = obj.neurites(idx,1).somad - d;
%                         obj.neurites(idx,1).crow = [obj.neurites(idx,1).crow, fRmax];%append new rownumber
%                         obj.neurites(idx,1).branch = [obj.neurites(idx,1).branch, CSVfile.Order(fR)]; %append new branch
%                     end
%                 end
               % idx = idx+1;
%             end
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


