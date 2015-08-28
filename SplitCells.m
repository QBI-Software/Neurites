function [cell1,cell2] = SplitCells(I,bgcolor)
% Function that performs color image thresholding using
% histogram derived. 
% Inputs:
% input - Input image as RGB
% bgcolor – value for bg color eg 256 for white, 0 for black.
% Output:
% 2 binary images of cells: cell1,cell2
    IG = rgb2gray(I);
       
    %histogram 
    [px,levels] = imhist(IG);
    figure
    hold on;
    bar(px);
    %histogram - exclude bg ?do nothing if not w or b
    min = 0;
    max = 256;
    if (bgcolor < 10)
        min = bgcolor + 1;
    elseif(bgcolor > 250)
        max = bgcolor - 1;
    end
    xlim([min max]);
    xlabel('color values')
    ylabel('pixels')
    grid on;
    
    title('Image Histogram');
   
    colors = [];
    i = 1;
    for pos = 1:254
        if px(pos) > 500
            colors(i,:) = [pos-1 px(pos)]
            i = i+1
        end
    end
    colors
    %Determine peaks from hist
    if (length(colors) > 0)
        c1 = sortrows(colors,2);
        peak1 = c1(end,:)
        peak2 = c1(end-1,:)
        peak3 = c1(end-2,:)
        
        cell1 = roicolor(IG,peak1(1));
        cell2 = roicolor(IG,peak2(1));
        othercolor = roicolor(IG, peak3(1));
        %Popup figure - can't get it to all display in box
        
        subplot(2,2,1)
        imshow(cell1)
        title('Cell 1')
        subplot(2,2,2)
        imshow(cell2)
        title('Cell 2')
        subplot(2,2,3)
        imshow(othercolor)
        title('Other')
    else
        error('Error: Unable to create histogram');
    end
    hold off;
end