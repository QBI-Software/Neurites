function fR=findCSVIndex(xC,yC,StartX,StartY,tol)
% Returns CSV Table row index matching xy coords    
    d = 0;
    fR = find((StartX >= xC-tol) & (StartX <=xC+tol));
    %if only one val
    if (isempty(fR))
        fR = find((StartY >= yC-tol) & (StartY <=yC+tol));
    end
    
    if (~isempty(fR) && (length(fR) > 1))
        P = [xC,yC]; 
        T = [StartX(fR) StartY(fR)];
        u=tol+1;
        s=0;
        for i=1:length(T)
            X = cat(1,P,T(i,:));
            d = pdist(X,'euclidean');
            if (d <= tol && u > d)
                    u = d;
                    s = fR(i);
            end
        end
        if (s > 0)
            fR =s;
        end
    end
    %If still empty or multiple
     if (isempty(fR) || length(fR) > 1)
         fR = [];
     end
    
end

