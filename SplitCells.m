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
    minm = min(levels);%0;
    maxm = max(levels);%256;
    if (bgcolor ==minm)
        minm = bgcolor + 1;
    elseif(bgcolor ==maxm)
        maxm = bgcolor - 1;
    end
    xlim([minm maxm]);
    xlabel('color values')
    ylabel('pixels')
    grid on;
    hold on;
    title('Image Histogram');
    
    colors = [];
    i = 1;
    for pos = 1:254
        if px(pos) > 200
            colors(i,:) = [pos-1 px(pos)];
            text(pos -0.4, px(pos) + 0.4, num2str(pos), 'VerticalAlignment', 'top', 'FontSize', 8)
            i = i+1;
        end
    end
    hold off;
    %colors
    %Determine peaks from hist
    if (length(colors) >= 2)
        c1 = sortrows(colors,2);
        peak1 = c1(end,:);
        peak2 = c1(end-1,:);
        
        %allow range
        p1 = [peak1(1)-10 : peak1(1)+10];
        p2 = [peak2(1)-10 : peak2(1)+10];
        cell1 = roicolor(IG,p1);
        cell2 = roicolor(IG,p2);
        subplot(2,2,1)
        imshow(cell1)
        title('Cell 1')
        subplot(2,2,2)
        imshow(cell2)
        title('Cell 2')
        if (length(colors) == 3)
            peak3 = c1(end-2,:);
            p3 = [peak3(1)-10 : peak3(1)+10];
            othercolor = roicolor(IG, p3);
            subplot(2,2,3)
            imshow(othercolor)
            title('Other')
        end
    else
        error('Error: Unable to create histogram');
    end
    hold off;
end