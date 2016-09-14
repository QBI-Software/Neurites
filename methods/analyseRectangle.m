function [annulus_area,neurites_area,regionMap] = analyseRectangle(Scale, Shiftx,Shifty,mask, maskedImage,centroid, direction,single,CSVFile )
    % Whole annulus, then partition via set arc length at different angles
    % Get masked image
    imshow(maskedImage); 
    hold on;
    scale = Scale * Scale;
    % Calculate full rectangle area
    mask_area = bwarea(mask)
    annulus_area= bwarea(maskedImage) %Total On pixels
    neurites_area = (mask_area - annulus_area)/scale %neurites only
    %remove annulus
    im2 = mask - maskedImage;
    %get params from mask
    stats = regionprops(mask,'BoundingBox');
    x = stats.BoundingBox(1);
    y = stats.BoundingBox(2);
    rheight = stats.BoundingBox(4);
    [maxx,maxy] = size(maskedImage);
    %Provide segmentation
    %arclength=45; %degrees
    if single
        a = y; 
    else
        %Ycoords: y = linspace(miny,maxy,by_y);
        %TODO: Determine maxy and provide option for direction: -x x -y y
        %from centroid
        switch direction
            case 'up'
                a = linspace(0,y,y/rheight); %upward
            case 'down'
                a = linspace(y,maxy,(maxy-y)/rheight); %downward
            case 'left'
                a = linspace(0,x,x/rwidth); %left
            case 'right'
                a = linspace(x,maxx,(maxx-x)/rwidth); %right
        end
    end
        %a = fliplr(a);
    %initialise
    stimregions = {length(a)};

    % Loop through segments
    %[h, w] = size(mask);%12390,13000
    colors=['b' 'g' 'r' 'c' 'm' 'y'];

    for i = 1:length(a)        
        m1 = mask; %copy
        im1 = im2;
        %Arc(i+1) = a(i);
        %Calculate mask for segment
        idx = mod(i,length(colors)-1);
        color = colors(idx+1);
        stats.BoundingBox(2)=a(i); %increment in y-axis
        rectangle('Position',stats.BoundingBox, 'EdgeColor',colors(idx+1))
        %P = plotArc(a(i),arclength, centroid(1),centroid(2),radius);
        xdata = linspace(stats.BoundingBox(1),stats.BoundingBox(1)+stats.BoundingBox(3));
        ydata = linspace(stats.BoundingBox(2),stats.BoundingBox(2)+stats.BoundingBox(4));
        %p1 = poly2mask(P.XData,P.YData,h,w);
        p1 = poly2mask(xdata,ydata,maxx,maxy);
        im1(~p1)=0;
        m1(~p1)=0;
        A= bwarea(m1)/scale;
        NA= bwarea(im1)/scale;
        imshow(im1)
        [B,~] = bwboundaries(im1,4,'noholes');
        
        N = countNeurites(im1,1,color);
        %show region boundary
        no = countNeurites(m1,1,color);
        
        %Identify neurites and Save data
        n = NeuritesStimulusRegion(a(i), A, i, NA, B, color);
        n = n.analyseBoundaries(Scale, Shiftx,Shifty,CSVFile);
        stimregions{i} = n;
       % waitbar(i/steps);
    end
   % close(hwb);
   %Create lookup map
   regionMap = containers.Map(a,stimregions);
    
end

