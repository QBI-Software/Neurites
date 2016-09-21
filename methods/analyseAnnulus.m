function [status,regionMap] = analyseAnnulus(Scale, Shiftx,Shifty,N, outer_radius,inner_radius,midline,arclength,single,CSVFile )
    % Whole annulus, then partition via set arc length at different angles
    %Load analyser: N.Iroi, N.maskedI, N.soma.centroid,N = N.loadMask(I);
    mask = N.Iroi;
    maskedImage = N.maskedI;
    centroid = N.soma.centroid;
    %Load DigiNeuron 
    hold on;
    scale = Scale * Scale;
    neuron = loadCSVTree(CSVFile,N.soma,Scale);
    status = sprintf('Loaded neuron with %d trees', length(neuron));
    disp(status)
    % Calculate full annulus area
    mask_area = bwarea(mask);
    annulus_area= bwarea(maskedImage); %Total On pixels
    neurites_area = (mask_area - annulus_area); %neurites only
    status = sprintf('Annulus with total %0.2f um2 containing neurites of %0.2f um2 (%0.5f percentage)', ...
        annulus_area/scale, neurites_area/scale, neurites_area/annulus_area);
    disp(status)
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
        P = plotArc(a(i),arclength, centroid(1),centroid(2),outer_radius);
        p1 = poly2mask(P.XData,P.YData,h,w);
        
        im1(~p1)=0;
        m1(~p1)=0;
        A= bwarea(m1)/scale;
        NA= bwarea(im1)/scale;
        imshow(im1)
        
        idx = mod(i,length(colors)-1);
        color = colors(idx+1);
        %countNeurites(im1,1,color);
        %show region boundary
        %countNeurites(m1,1,color);
        % determine boundaries
        %[B,L,N,A] = bwboundaries(m1,4,'noholes');
        %Border points
        [B,~] = bwboundaries(m1,4,'noholes'); %mask outline
        b = B{1};
        plot(b(:,2),b(:,1),color,'LineWidth',2)
        xdata= b(:,2);
        ydata = b(:,1);
        pary = [xdata ydata];  
        
        %Identify neurites and Save data
        [dcache,FR] = findCSVdata(CSVFile,Scale,Shiftx,Shifty,pary,xdata,ydata,color);
        n = NeuritesStimulusRegion(i, a(i), A, arclength, NA, im1, color,'annulus',neuron);
        %n = NeuritesStimulusRegion(i,a(i), A, arclength, NA, B, color, 'annulus');
        %n = n.analyseBoundaries(Scale, Shiftx,Shifty,CSVFile);
        n = n.setScale(Scale, Shiftx,Shifty);
        n = n.analyseROI(CSVFile,FR,dcache);
        stimregions{i} = n;
       % waitbar(i/steps);
    end
   % close(hwb);
   %Create lookup map
   regionMap = containers.Map(a,stimregions);
    
end

