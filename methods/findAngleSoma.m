function [theta,rho,deg] = findAngleSoma(x,y,soma,scale)
%relative to somacentroid 
    dx = x-soma.centroid(:,1); %x1-x0
    dy = y-soma.centroid(:,2); %y1-y0
    theta= atan2(dy,dx); %radians flipped by 180o == pi
    rho = sqrt(dx.^2 + dy.^2)/scale;
    deg = radtodeg(theta); 
    if (deg < 0)
        deg = deg + 360; %negatives
    end
end
