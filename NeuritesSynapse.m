classdef NeuritesSynapse
   properties
      
      %Binary img analysis
      SegmentC1 %cell 1 neurite
      SegmentC2 %cell2 neurite
      RegionC1 %synaptic region relative to cell 1
      RegionC2 %synaptic region relative to cell 2
      MedianC1
      MedianC2
      SynapseType %1=enpassant 2=end terminal 3=intersection
      %Mapping to csv
      shiftx  % x-shift to world co-ords
      shifty  % y-shift to world co-ords
      scale
      StartXY1 % Start coords from csv for cell1
      StartXY2 % Start coords from csv for cell2
      %neuron measurements
      NeuriteLengthC1
      NeuriteLengthC2
      NeuriteEndC1
      NeuriteEndC2
      DistanceC1
      DistanceC2
      BranchPointC1
      BranchPointC2
      BranchTypeC1
      BranchTypeC2
      EndpointsC1
      EndpointsC2
      SomaC1
      SomaC2
      SomapointsC1
      SomapointsC2
      fR1
      fR2
      
   end
   methods
      
      function obj = NeuritesSynapse(redSegment, greenSegment, redhighlight, greenhighlight, ...
           rmarkerx,rmarkery,gmarkerx,gmarkery,type)
        obj.SegmentC1 = redSegment;
        obj.SegmentC2 = greenSegment;
        obj.RegionC1 = redhighlight;
        obj.RegionC2 = greenhighlight;
        obj.MedianC1 = [rmarkerx,rmarkery];
        obj.MedianC2 = [gmarkerx,gmarkery];
        obj.SynapseType = type;
      end
      function [x,y] = img2Coords(obj,imgx, imgy)
        x = (imgx - obj.shiftx)/obj.scale;
        y = (-imgy - obj.shifty)/obj.scale;

      end
      function [imgx,imgy] = coords2Img(obj,x, y)
        imgx = (x * obj.scale) + obj.shiftx;
        imgy = -((y * obj.scale) + obj.shifty);

      end
      %isequal based on measurements only
      function rslt = isequal(obj,synA)
          rslt = 0;
          %isa(x,'className');
          if (isa(synA,'NeuritesSynapse'))
              if (~isempty(obj.NeuriteLengthC1))
                if ((synA.NeuriteLengthC1 == obj.NeuriteLengthC1) && ...
                  (synA.NeuriteLengthC2 == obj.NeuriteLengthC2) && ... 
                  (synA.DistanceC1 == obj.DistanceC1) && ...
                  (synA.DistanceC2 == obj.DistanceC2) && ...
                  (synA.BranchPointC1 == obj.BranchPointC1) && ...
                  (synA.BranchPointC2 == obj.BranchPointC2) )
                  rslt = 1;
                end
              end
          end
      end
      %remove a branch containing point xy
      function obj = removeBranch(obj,c,T,x,y)
          sx = (x - obj.shiftx)/obj.scale;
          sy = (-y - obj.shifty)/obj.scale;
          xi = find(T.StartX == sx);
          yi = find(T.StartY == sy);
          ci = intersect(xi,yi,'rows');
          xpoints=[];
          ypoints=[];
          if (c==1)
              somapoints = obj.SomapointsC1;
              soma = obj.SomaC1;
          else
              somapoints = obj.SomapointsC2;
              soma = obj.SomaC2;
          end
          for (i=1:length(ci)) 
              fR = ci(i);
              tree = T.Tree(fR);
              order = T.Order(fR);
              rows = find(T.Tree == tree & T.Order == order);
              for(j=1:length(rows))
                  xpoints(end+1) = (T.StartX(rows(j)) * obj.scale) + obj.shiftx;
                  ypoints(end+1) = -((T.StartY(rows(j)) * obj.scale) + obj.shifty);
              end
          end
          %remove
          if (length(xpoints) > 0)
              saved = [];
              removed = [];
              for (j = 1:length(somapoints))
                if (any(bsxfun(@eq,xpoints(:),somapoints(j,1))) && any(bsxfun(@eq,ypoints(:),somapoints(j,2))))
                    removed(end+1) = j;
                    %somapoints(j,:) = [];
                else
                    saved = cat(1,saved,somapoints(j,:));
                end
              end
              somapoints = saved;
              %recalculate synapse
              soma = obj.calculateDistance(somapoints);
              %save back
              if (c==1)
                obj.SomapointsC1 = somapoints;
                obj.SomaC1 = soma;
              else
                obj.SomapointsC2 = somapoints;
                obj.SomaC2 = soma;
              end
          end
      end %end function
      
      function dist = calculateDistance(obj,somapoints)
          d = 0;
          for i=1:length(somapoints)-1
              pts = cat(1, somapoints(i,:), somapoints(i+1,:));
              d = d + pdist(pts,'euclidean');
          end
          dist = d/obj.scale;    
      end
   end
   
end

