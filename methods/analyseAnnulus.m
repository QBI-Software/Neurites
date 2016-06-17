function [T,colnames] = analyseAnnulus(Scale, mask, maskedImage, centroid, radius,midline,arclength,show,single )
    % Whole annulus, then partition via set arc length at different angles
    % Get masked image
    imshow(maskedImage); 
    hold on;
    scale = Scale * Scale;
    % Calculate full annulus area
    mask_area = bwarea(mask)
    annulus_area= bwarea(maskedImage) %Total On pixels
    neurites_area = (mask_area - annulus_area)/scale %neurites only
    %remove annulus
    im2 = mask - maskedImage;
    %Provide segmentation
    %arclength=45; %degrees
    if single
        a = midline; 
    else
        a = linspace(arclength,360,360/arclength); %[360;45;90;135;180;225;270;315];
    end
        %a = fliplr(a);
    %initialise
    %Arc = zeros(length(a)+1);
    %A = zeros(length(a)+1);
    %NA = zeros(length(a)+1);
    %N = zeros(length(a)+1);
    %First value is whole annulus
    Arc(1) = 0;
    A(1)=mask_area/scale;
    NA(1)=neurites_area;
    N(1) = countNeurites(im2,0,'w');
    C(1) = 0;
    BBX(1) = 0;
    BBY(1) = 0;
    BBW(1) = 0;
    BBH(1)= 0;
    BBArea(1) = 0;
    % Loop through segments
    [h, w] = size(mask);%12390,13000
    colors=['b' 'g' 'r' 'c' 'm' 'y'];
   % hwb = waitbar(0,'Running analysis...','WindowStyle','modal');
    steps = length(a);
        
    for i = 1:length(a)
        
        m1 = mask; %copy
        im1 = im2;
        Arc(i+1) = a(i);
        %Calculate mask for segment
        P = plotArc(a(i),arclength, centroid(1),centroid(2),radius);
        p1 = poly2mask(P.XData,P.YData,h,w);
        im1(~p1)=0;
        m1(~p1)=0;
        A(i+1)= bwarea(m1)/scale;
        NA(i+1)= bwarea(im1)/scale;
        imshow(im1)
        idx = mod(i,length(colors)-1);
        color = colors(idx+1);
        N(i+1) = countNeurites(im1,1,color);
        C(i+1) = idx+1;
        %show region boundary
        no = countNeurites(m1,1,color);
        st = regionprops(m1, 'BoundingBox' );
        BBX(i+1) = st.BoundingBox(1);
        BBY(i+1) = st.BoundingBox(2);
        BBW(i+1) = st.BoundingBox(3)/scale;
        BBH(i+1) = st.BoundingBox(4)/scale;
        BBArea(i+1) = (st.BoundingBox(3) * st.BoundingBox(4))/scale;
       % waitbar(i/steps);
    end
   % close(hwb);
    ArcMidline = Arc';
    Area = A';
    NeuritesArea = NA';
    Number = N';
    Color = C';
    BBx = BBX';
    BBy = BBY';
    BBw = BBW';
    BBh = BBH';
    BBarea = BBArea';
    T = table(ArcMidline,Area,NeuritesArea,Number,Color,BBx,BBy,BBw,BBh,BBarea)
    colnames = {'ArcMidline' 'Area' 'NeuritesArea' 'Number' 'Color' ...
        'BBx' 'BBy' 'BBw' 'BBh' 'BBarea'};
    %Show data
    if show
        htable = findobj('Tag','uitableResults');
        set(htable,'data',[T.ArcMidline T.Area T.NeuritesArea T.Number, T.Color ...
            T.BBx T.BBy T.BBw T.BBh T.BBarea],'ColumnName', colnames);
    end
    
end

