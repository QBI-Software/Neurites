function [output] = RGBThreshold(input,thresholds)
% Function that performs color image thresholding using
% user-defined threshold values. Emphasises red areas.
% Inputs:
% input - Input image
% thresholds – 1x3 matrix with the threshold values
% for the R, G and B color channels.
% Output:
% output - Output image (binary)
    red_bin = input(:,:,1) > thresholds(1); % Red thresholding
    green_bin = input(:,:,2) < thresholds(2); % Green thresholding
    blue_bin = input(:,:,3) < thresholds(3); % Blue thresholding
    output = red_bin & green_bin & blue_bin; % Final image
end