classdef NeuritesSoma
    properties
        area
        perimeter
        centroid
    end
    methods
        function obj = NeuritesSoma(area,perimeter,centroid)
            obj.area = area;
            obj.perimeter = perimeter;
            obj.centroid = centroid;       
          end
    end
    
end