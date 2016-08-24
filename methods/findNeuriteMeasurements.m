%Find cumulative length, volume, sa of neurite to end or to finalxy point
%struct: 
function [d,v,s,endpoints, branchpoints] = findNeuriteMeasurements(syn,x,y,finalxy,T,fR)
    maxfr = height(T);
    endxy = [T.EndX(fR), T.EndY(fR)];
    synxy = [x,y];
    xpoints = [(x * syn.scale) + syn.shiftx];
    ypoints = [-((y  * syn.scale) + syn.shifty)];
    d = distanceBetween(synxy,endxy,'euclidean');
    radius = T.AvgDiameter__m_(fR)/2;
    v = pi * d * radius^2; 
    s = (2 * pi * radius * d) + (2 * pi * radius^2); 
    cont = 1;
    branchpoints = [fR];
    while(strcmp(T.PointType(fR),'EP')== 0 && fR < maxfr && cont)
        if (strcmp(T.PointType(fR),'BP')> 0)
            branchpoints(end+1)=fR;
        end
        d = d + T.Length__m_(fR);
        v = v + T.Volume__m__(fR);
        s = s + T.SurfaceArea__m__(fR);
        %Save points for display
        endx = (T.StartX(fR) * syn.scale) + syn.shiftx
        endy = -((T.StartY(fR) * syn.scale) + syn.shifty)
        
        if (finalxy(1) > endx && finalxy(2) > endy)
            xpoints(end+1)=endx ;
            ypoints(end+1)=endy ;
        else
            cont = 0;
        end
        fR = fR + 1;
    end    
    endpoints =cat(2,xpoints',ypoints');
end

