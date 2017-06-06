%line1
x1  = [7.8 8.5];
y1  = [0.96 0.94];
%line2
x2 = [8.25 8.25];
y2 = [0 0.99];
%fit linear polynomial
p1 = polyfit(x1,y1,1);
p2 = polyfit(x2,y2,1);
%calculate intersection
x_intersect = fzero(@(x) polyval(p1-p2,x),3);
y_intersect = polyval(p1,x_intersect);
line(x1,y1);
hold on;
line(x2,y2);
plot(x_intersect,y_intersect,'r*')