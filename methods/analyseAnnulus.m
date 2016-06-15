function [T,colnames] = analyseAnnulus(Scale, mask, maskedImage, centroid, radius,show)
    % Whole annulus, then partition via set arc length at different angles
    % Get masked image
    imshow(maskedImage); 
    hold on;
    scale = Scale * Scale;
    % Calculate full annulus area
    mask_area = bwarea(mask)
    annulus_area= bwarea(maskedImage) %Total On pixels
    neurites_area = (mask_area - annulus_area)/scale %neurites only
    
    %Provide segmentation
    neurites_num= 5; %replace with method to calculate number neurites
    arclength=45; %degrees
    a = linspace(45,360,8); %[360;45;90;135;180;225;270;315];
    a = fliplr(a);
    %First value is whole annulus
    Arc(1) = 0;
    A(1)=mask_area/scale;
    NA(1)=neurites_area;
    N(1) = countNeurites(maskedImage,0);
    % Loop through segments
    [h w] = size(mask);%12390,13000
    colors=['b' 'g' 'r' 'c' 'm' 'y'];
    for i = 1:length(a)
        m1 = mask; %copy
        im1 = maskedImage;
        Arc(i+1) = a(i);
        %Calculate mask for segment
        P = plotArc(a(i),arclength, centroid(1),centroid(2),radius);
        p1 = poly2mask(P.XData,P.YData,h,w);
        im1(~p1)=0;
        m1(~p1)=0;
        A(i+1)= bwarea(m1)/scale;%polyarea(P.XData,P.YData)/Scale;
        NA(i+1)= (bwarea(m1) - bwarea(im1))/scale;
        imshow(im1)
        idx = mod(i,length(colors));
        if idx == 0
           idx = length(colors);
        end
        color = colors(idx);
        N(i+1) = countNeurites(im1,1,color);
        
    end
    ArcMidline = Arc';
    Area = A';
    NeuritesArea = NA';
    Number = N';    
    T = table(ArcMidline,Area,NeuritesArea,Number)
    colnames = {'ArcMidline' 'Area' 'NeuritesArea' 'Number'};
    %Show data
    if show
        htable = findobj('Tag','uitableResults');
        set(htable,'data',[T.ArcMidline T.Area T.NeuritesArea T.Number],'ColumnName', colnames);
    end
end

