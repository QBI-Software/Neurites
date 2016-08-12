function [annulus_area,neurites_area,regionMap] = analyseAnnulus(Scale, Shiftx,Shifty,mask, maskedImage, centroid, radius,midline,arclength,single,CSVFile )
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
    stimregions = {length(a)};

    % Loop through segments
    [h, w] = size(mask);%12390,13000
    colors=['b' 'g' 'r' 'c' 'm' 'y'];

    for i = 1:length(a)
        
        m1 = mask; %copy
        im1 = im2;
        %Arc(i+1) = a(i);
        %Calculate mask for segment
        P = plotArc(a(i),arclength, centroid(1),centroid(2),radius);
        p1 = poly2mask(P.XData,P.YData,h,w);
        im1(~p1)=0;
        m1(~p1)=0;
        A= bwarea(m1)/scale;
        NA= bwarea(im1)/scale;
        imshow(im1)
        [B,~] = bwboundaries(im1,4,'noholes');
        idx = mod(i,length(colors)-1);
        color = colors(idx+1);
        N = countNeurites(im1,1,color);
        %show region boundary
        no = countNeurites(m1,1,color);
        
        %Identify neurites and Save data
        n = NeuritesStimulusRegion(a(i), A, arclength, NA, B, color);
        n = n.analyseBoundaries(Scale, Shiftx,Shifty,CSVFile);
        stimregions{i} = n;
       % waitbar(i/steps);
    end
   % close(hwb);
   %Create lookup map
   regionMap = containers.Map(a,stimregions);
    
end

