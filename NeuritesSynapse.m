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
      DistanceC1
      DistanceC2
      BranchPointC1
      BranchPointC2
      BranchTypeC1
      BranchTypeC2
      
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
      %function to lookup neuron measurements from csv
   end
   
end

