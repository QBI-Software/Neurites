function [imgx,imgy] = coords2Img(scale,shiftx,shifty,x, y) 
%Converts CSV coords to IMAGE coords
%   
    imgx = (x * scale) + shiftx;
    imgy = -((y * scale) + shifty);

end

