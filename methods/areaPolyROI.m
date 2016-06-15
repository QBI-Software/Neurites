function [RCols, RData] = areaPolyROI( dataroi, scale )
% Calculate PolyROI area
%  returns a table columns and data

    %ROI stats
    numPixels = sum(dataroi.BW(:))
    roi_area = polyarea(dataroi.xi,dataroi.yi);
    
    %Convert ROI to um2 if csv config loaded
    roi_um2 = roi_area / scale;
    %Results table
    RCols={'ROI Number Pixels', 'ROI Area (px)','ROI Area (um2)'};
    RData=[numPixels, roi_area,roi_um2];
    

end

