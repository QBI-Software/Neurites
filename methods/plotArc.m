function P = plotArc(midline,length,cx,cy,r)
% Plot a circular arc as a pie wedge.
% midline in degrees - middle of arc
% length in degrees - full length of arc
% a is start of arc in radians, 
% b is end of arc in radians, 
% (cx,cy) is the center of the circle.
% r is the radius.
% Try this:   plot_arc(pi/4,3*pi/4,9,-4,3)
% Based on Author:  Matt Fig
m = deg2rad(midline);
l = deg2rad(length);
a = m - l/2;
b = m + l/2;
t = linspace(a,b);
x = r*cos(t) + cx;
y = r*sin(t) + cy;
x = [x cx x(1)];
y = [y cy y(1)];
P = fill(x,y,'g');
axis([cx-r-1 cx+r+1 cy-r-1 cy+r+1]) 
axis square;
if ~nargout
    clear P
end

