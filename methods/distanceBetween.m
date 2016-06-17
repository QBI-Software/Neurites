function d = distanceBetween(A, B, type)
    d= 0;
    if (isempty(type))
        type = 'euclidean';
    end
    X = cat(1,A,B);
    d = pdist(X,type);
end


