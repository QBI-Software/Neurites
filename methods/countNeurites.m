function number = countNeurites(im,showplots,color)
% Count neurites within a binary image
%   Uses boundaries
    [B,L] = bwboundaries(im,'noholes');
    %imshow(label2rgb(L, @jet, [.5 .5 .5]))
    %hold on
    
    if (showplots)
        for k = 1:length(B)
           boundary = B{k};
           plot(boundary(:,2), boundary(:,1), color, 'LineWidth', 2)
        end
    end
    number = length(B);
end

