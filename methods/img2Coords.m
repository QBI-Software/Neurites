function [x,y] = img2Coords(scale,shiftx,shifty,imgx, imgy)
    x = (imgx - shiftx)/scale;
    y = (-imgy - shifty)/scale;

end
