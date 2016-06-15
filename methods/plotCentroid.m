function p1 = plotCentroid( CentroidX, CentroidY )
%Plot Centroid on current fig
%   Detailed explanation goes here
gca
hold on
p1 = plot(CentroidX, CentroidY,'color', 'b',...
                'marker','x','linestyle','none','LineWidth', 2);
hold off
end

