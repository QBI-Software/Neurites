function [roi_area,neurites_area,regionMap] = analyseRectangle(Scale, Shiftx,Shifty,N, direction,single,CSVFile )
    %Load analyser: N.Iroi, N.maskedI, N.soma.centroid,N = N.loadMask(I);
    mask = N.Iroi;
    maskedImage = N.maskedI;
    %Load DigiNeuron
    neuron = loadCSVTree(CSVFile,N.soma,Scale);
    sprintf('Loaded neuron with %d trees', length(neuron));
    % Get masked image
    imshow(maskedImage); 
    hold on;
    scale = Scale * Scale;
    % Calculate full rectangle area
    mask_area = bwarea(mask)
    roi_area= bwarea(maskedImage) %Total On pixels
    neurites_area = (mask_area - roi_area) %neurites only
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
                L = rheight/Scale;
            case {'left', 'right'}
                stats.BoundingBox(1) = a(i);
                L = rwidth/Scale;
        end
        %Calculate mask for segment
        idx = mod(i,length(colors)-1);
        color = colors(idx+1);
        rectangle('Position',stats.BoundingBox, 'EdgeColor',color)
       
        %create ROI
        x0 = stats.BoundingBox(1);
        y0 = stats.BoundingBox(2);
        xdata = [x0,x0+rwidth,x0+rwidth,x0,x0]; %row vertices
        ydata =[y0,y0,y0+rheight,y0+rheight,y0]; %col vertices
        p1 = poly2mask(xdata,ydata,maxx,maxy);
        N = N.loadMask(p1);
        im1 = N.maskedI;
        m1 = N.Iroi;
        im1 = m1 - im1; %invert
        A= bwarea(m1)/scale;
        NA= bwarea(im1)/scale;
        imshow(im1)

        %Border points
        asize = 1000;
        px = linspace(x,x+rwidth,asize);
        py = linspace(y,y+rheight,asize);
        pleft = [linspace(x,x,asize)',py']; 
        pright = [linspace(x+rwidth,x+rwidth,asize)',py'];
        ptop = [px',linspace(y,y,asize)'];
        pbottom = [px',linspace(y+rheight,y+rheight,asize)'];
        
        %Show borders
        plot(pleft(:,1), pleft(:,2),color)
        plot(pright(:,1), pright(:,2),color)
        plot(ptop(:,1), ptop(:,2),color)
        plot(pbottom(:,1), pbottom(:,2),color)
        pary = [pleft pright ptop pbottom];      

        %Find and analyse CSV data
        [dcache,FR] = findCSVdata(CSVFile,Scale,Shiftx,Shifty,pary,xdata,ydata,color)
        %Identify neurites and Save data (id, midline, sarea, slength, neuritesarea, boundaries, color, type)
        n = NeuritesStimulusRegion(i, a(i), A, L, NA, im1, color,'rect',neuron);
        n = n.setScale(Scale, Shiftx,Shifty);
        n = n.analyseROI(CSVFile,FR,dcache);
        stimregions{i} = n;
       % waitbar(i/steps);
    end
   % close(hwb);
   %Create lookup map
   regionMap = containers.Map(a,stimregions);
    
end

