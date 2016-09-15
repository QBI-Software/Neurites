function [roi_area,neurites_area,regionMap] = analyseRectangle(Scale, Shiftx,Shifty,N, direction,single,CSVFile )
    %Load analyser: N.Iroi, N.maskedI, N.soma.centroid,N = N.loadMask(I);
    mask = N.Iroi;
    maskedImage = N.maskedI;
    
    
    % Get masked image
    imshow(maskedImage); 
    hold on;
    scale = Scale * Scale;
    % Calculate full rectangle area
    mask_area = bwarea(mask)
    roi_area= bwarea(maskedImage) %Total On pixels
    neurites_area = (mask_area - roi_area)/scale %neurites only
    %remove annulus
    %im2 = mask - maskedImage;
    %get params from mask
    stats = regionprops(mask,'BoundingBox');
    x = stats.BoundingBox(1);
    y = stats.BoundingBox(2);
    rwidth = stats.BoundingBox(3);
    rheight = stats.BoundingBox(4);
    [maxx,maxy] = size(maskedImage);
    
    if single
        a = y; 
    else
        a = [];
        switch direction
%             case 'up'
%                 a = []linspace(y,rheight,(y+rheight)/rheight); %upward stats.BoundingBox(2)-rheight
%             case 'down'
%                 a = linspace(y,maxy,(maxy-y+rheight)/rheight); %downward stats.BoundingBox(2)+rheight
%             case 'left'
%                 a = linspace(x,0,(x+rwidth)/rwidth); %left stats.BoundingBox(1)-rwidth
%             case 'right'
%                 a = linspace(x,maxx-rwidth,(maxx-x+rwidth)/rwidth); %right stats.BoundingBox(1)+rwidth
            case 'up'
                for j=1:(y/rheight)
                    a(j) = y;
                    y = y-rheight; %upward stats.BoundingBox(2)-rheight
                end
            case 'down'
                for j=1:((maxy-y)/rheight)
                    a(j) = y;
                    y = y+rheight; 
                end
            case 'left'
                for j=1:(x/rwidth)
                    a(j) = x;
                    x = x-rwidth; %upward stats.BoundingBox(2)-rheight
                end
            case 'right'
                for j=1:((maxx-x)/rwidth)
                    a(j) = x;
                    x = x+rwidth; 
                end

        end
    end
        %a = fliplr(a);
    %initialise
    stimregions = {length(a)};

    % Loop through segments
    %[h, w] = size(mask);%12390,13000
    colors=['b' 'g' 'r' 'c' 'm' 'y'];

    for i = 1:length(a)        
        %increment
        switch direction
            case {'up', 'down'}
                stats.BoundingBox(2) = a(i);%stats.BoundingBox(2)-rheight
                L = rheight;
            case {'left', 'right'}
                stats.BoundingBox(1) = a(i);
                L = rwidth;
        end
        %Calculate mask for segment
        idx = mod(i,length(colors)-1);
        color = colors(idx+1);
        rectangle('Position',stats.BoundingBox, 'EdgeColor',colors(idx+1))
       % hr = imrect(gca,stats.BoundingBox)
        %create ROI
        x0 = stats.BoundingBox(1);
        y0 = stats.BoundingBox(2);
        xdata = [x0,x0+rwidth,x0+rwidth,x0,x0]; %row vertices
        ydata =[y0,y0,y0+rheight,y0+rheight,y0]; %col vertices
        p1 = poly2mask(xdata,ydata,maxx,maxy);
        N = N.loadMask(p1);
        im1 = N.maskedI;
        m1 = N.Iroi;
        im1 = m1 - im1; %inverse
        A= bwarea(m1)/scale;
        NA= bwarea(im1)/scale;
        imshow(im1)
        [B,~] = bwboundaries(im1,4,'noholes');
        %show neurites
        countNeurites(im1,1,color);
        %show region boundary
        countNeurites(m1,1,color);
        
        %Identify neurites and Save data (id, midline, sarea, slength, neuritesarea, boundaries, color, type)
        n = NeuritesStimulusRegion(i, a(i), A, L, NA, B, color,'rect');
        n = n.analyseBoundaries(Scale, Shiftx,Shifty,CSVFile);
        stimregions{i} = n;
       % waitbar(i/steps);
    end
   % close(hwb);
   %Create lookup map
   regionMap = containers.Map(a,stimregions);
    
end

