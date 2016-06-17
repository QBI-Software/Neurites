function [d,endpoints] = findLength2NeuriteEnd(syn,x,y,T,fR)
    maxfr = height(T);
    endxy = [T.EndX(fR), T.EndY(fR)];
    synxy = [x,y];
    xpoints = [(x * syn.scale) + syn.shiftx];
    ypoints = [-((y  * syn.scale) + syn.shifty)];
    d = distanceBetween(synxy,endxy,'euclidean');
    while(strcmp(T.PointType(fR),'EP')== 0 && fR < maxfr)
        d = d + T.Length__m_(fR);
        %Save points for display
        xpoints(end+1) = (T.StartX(fR) * syn.scale) + syn.shiftx;
        ypoints(end+1) = -((T.StartY(fR) * syn.scale) + syn.shifty);
        fR = fR + 1;
    end    
    endpoints =cat(2,xpoints',ypoints');
end

