function d = findDistanceToSoma(x,y,startxy,lengthxy)
    d = distanceBetween([x,y],startxy,'euclidean');
    d = lengthxy - d;
end

