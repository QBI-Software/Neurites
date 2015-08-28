classdef NeuritesSynapse
   properties
      x  %world co-ords
      y
      bx %binary img x,y
      by
      %Binary img analysis
      SegmentC1 %cell 1 neurite
      SegmentC2 %cell2 neurite
      RegionC1 %synaptic region relative to cell 1
      RegionC2 %synaptic region relative to cell 2
      MedianC1
      MedianC2
      SynapseType %1=enpassant 2=end terminal 3=intersection
      %neuron measurements
      NeuriteLengthC1
      NeuriteLengthC2
      DistanceC1
      DistanceC2
      BranchPointC1
      BranchPointC2
      TotalBranchPtsC1
      TotalBranchPtsC2
      
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
      function [x,y] = getMedianWcoords(imgx, imgy)
        x = imgx;
        y = imgy;

      end
      %function to lookup neuron measurements from csv
   end
   
end

