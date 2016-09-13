function outermask = applyRectangle(w,h, scale, centroid)
        rw = (h * scale);
        rh = (w * scale);
        x = centroid(1) - (rw/2);  % min_x
        y = centroid(2) - (rh/2);  % min_y
        outer = [x y rw rh]; % [min_x min_y diam diam]
        % Plot center
        
        hold on;
        plot(centroid(1), centroid(2), 'r+', 'LineWidth', 2, 'MarkerSize', 5);
        grid on;
        %setPosition(h,pos) sets the rectangle h to a new position. The new position, pos, has the form [xmin ymin width height].
        hRect = imrect(gca, outer);
        status = sprintf('Move by dragging then Double-Click on rectangle to accept position.');
        msgbox(status);
        position = wait(hRect); %confirm position
        outermask = hRect.createMask;
       
end

