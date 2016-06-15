function [Scale, Shiftx, Shifty] = estimateOverlay(handles, I,cell1csv)
%Estimate parameters for overlay of CSV data onto image
    Img = rgb2gray(I);
    Tb = readtable(cell1csv);
    
    %CSV data in um
    mi = min(Tb.StartX);
    if (min(Tb.EndX) < mi)
        rowi = find(Tb.EndX == min(Tb.EndX));
        S1 = [Tb.EndX(rowi),Tb.EndY(rowi)];
    else
        rowi = find(Tb.StartX == min(Tb.StartX));
        S1 = [Tb.StartX(rowi),Tb.StartY(rowi)];
    end
    mx = max(Tb.StartX);
    if (max(Tb.EndX) > mx)
        rowi = find(Tb.EndX == max(Tb.EndX));
        S2 = [Tb.EndX(rowi),Tb.EndY(rowi)];
    else
        rowi = find(Tb.StartX == max(Tb.StartX));
        S2 = [Tb.StartX(rowi),Tb.StartY(rowi)];
        
    end
    if (length(S1(:,1)) > 1)
        rowi = find(S1(:,2) == max(S1(:,2)));
        S1 = [S1(rowi,1),S1(rowi,2)];
    end
    if (length(S2(:,1)) > 1)
        rowi = find(S2(:,2) == max(S2(:,2)));
        S2 = [S2(rowi,1),S2(rowi,2)];
    end
    C1 = cat(1,S1,S2);
    csvlength = pdist(C1,'euclidean');
    csvlength = csvlength + (0.1*csvlength); %adjust for slightly shorter traces
    %Image coords - 
    LI = bwlabel(Img);
    s = regionprops(LI, 'Area', 'BoundingBox','Extrema');
    i = 1;
    Area = cat(1,s.Area);
    if (length(s) > 1)
        i = find(Area==max(Area));
    end
%     S3=s(i).Extrema(1,:);
%     S4=s(i).Extrema(4,:);
%     S3=[s(i).BoundingBox(1) s(i).BoundingBox(2)];
%     S4=[s(i).BoundingBox(1)+s(i).BoundingBox(3)  s(i).BoundingBox(2)+s(i).BoundingBox(4)];
    Ic = corner(Img,'Harris','SensitivityFactor',0.02);
    xi = Ic(:,1)==min(Ic(:,1));
    S3 = Ic(xi,:)
    xm = Ic(:,1)==max(Ic(:,1));
    S4 = Ic(xm,:)
    if (length(S3(:,1)) > 1)
        rowi = find(S3(:,2) == max(S3(:,2)));
        S3 = [S3(rowi,1),S3(rowi,2)];
    end
    if (length(S4(:,1)) > 1)
        rowi = find(S4(:,2) == max(S4(:,2)));
        S4 = [S4(rowi,1),S4(rowi,2)];
    end
    % Find scale
    C2 = cat(1,S3,S4);
    imglength = pdist(C2,'euclidean');
    Scale = round(imglength/csvlength,1)
    S = S1 * Scale;
    
    %Plot
    axes(handles.axes1);
    imshow(Img);
    hold on;
    %plot(S(1),S(2),'color','b','Marker','o','LineWidth', 2);
    plot(S3(1),S3(2),'color','y','Marker','x','LineWidth', 2);
    plot(S4(1),S4(2),'color','r','Marker','x','LineWidth', 2);
    hold off;
    %Load data
    
    Shiftx = round(S3(1) - S(1),2);
    Shifty = round(-S3(2) - S(2),2);
    
end

