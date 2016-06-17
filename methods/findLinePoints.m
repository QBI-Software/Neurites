function [p1,p2] = findLinePoints(region)
% Finds end points on a line specified by region (r,c)
    rx = region(:,2);
    ry = region(:,1);
    [rxmin,ri1] = min(rx);
    [rxmax,ri2] = max(rx);
    [rymin,ri3] = min(ry);
    [rymax,ri4] = max(ry);
    if (rxmax - rxmin > 0)
        p1 = [rxmin, ry(ri1)];
        p2 = [rxmax, ry(ri2)];
    else
        p1 = [rx(ri3), rymin ];
        p2 = [rx(ri4), rymax ];
    end

end

