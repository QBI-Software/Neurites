function c = findLength(region, resolution)
%Find length of region - assumes straight line 
    if(isempty(region))
        c = 0;
    else
        sortedregion = sortrows(region,1);
        top = sortedregion(1,:);
        bottom = sortedregion(end,:);
        x1 = top(1);
        y1 = top(2);
        x2 = bottom(1);
        y2 = bottom(2);
        V = [x1 y1;x2 y2];
        c = pdist(V,'euclidean') * resolution;
    end
end


