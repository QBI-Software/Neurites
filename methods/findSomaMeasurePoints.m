%Get points to soma,and measurements: distance, volume, surf area
function [d,v,s,points] = findSomaMeasurePoints(syn,x,y,T,fR)
    maxfr = height(T);
    startxy = [T.StartX(fR), T.StartY(fR)];
    synxy = [x,y];
    tree = T.Tree(fR);
    order = T.Order(fR);
    cacheTree = tree;
    cacheOrder = order;
    %endpoints = {}
    xpoints = [];
    ypoints = [];
    d = distanceBetween(synxy,startxy,'euclidean');
    radius = T.AvgDiameter__m_(fR)/2;
    v = pi * d * radius^2; %0;
    s = (2 * pi * radius * d) + (2 * pi * radius^2); %0
    cacheStart = T.StartX(fR);
    while(tree == cacheTree  && fR < maxfr && fR > 0)
        %check tree and branch order strcmp(T.PointType(fR),'EP')== 0 && 
        tree = T.Tree(fR);
        order = T.Order(fR);
        if (cacheStart == T.EndX(fR) && (order == cacheOrder || order == cacheOrder-1))
            cacheStart = T.StartX(fR);
            d = d + T.Length__m_(fR);
            v = v + T.Volume__m__(fR);
            s = s + T.SurfaceArea__m__(fR);
            %Save points for display
            xpoints(end+1) = (T.StartX(fR) * syn.scale) + syn.shiftx;
            ypoints(end+1) = -((T.StartY(fR) * syn.scale) + syn.shifty);
        
           if (order == cacheOrder-1 )
                cacheOrder = order;
           end
        end
       fR = fR - 1;
    end    
    points =cat(2,xpoints',ypoints');
end
