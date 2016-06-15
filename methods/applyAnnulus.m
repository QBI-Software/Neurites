function outermask = applyAnnulus(id, od, scale, centroid)
        rid = (id * scale);
        rod = (od * scale);
        x = centroid(1) - (rod/2);  % min_x
        y = centroid(2) - (rod/2);  % min_y
        outer = [x y rod rod]; % [min_x min_y diam diam]
        % Plot center
        
        hold on;
        plot(centroid(1), centroid(2), 'r+', 'LineWidth', 2, 'MarkerSize', 5);
        grid on;
        hCircle1 = imellipse(gca, outer);
        status = sprintf('Double-Click on OUTER annulus to accept position.');
        msgbox(status);
        position = wait(hCircle1); %confirm position
        outermask = hCircle1.createMask;
        %[xmin ymin width height] = position
        %[3995.90774907749 4541.9667896679 3500 3500]
        x = centroid(1) - (rid/2)
        y = centroid(2) - (rid/2)
        inner = [x y rid rid]
        hCircle2 = imellipse(gca, inner);
        status = sprintf('Double-Click on INNER annulus to accept position.');
        msgbox(status);
        position = wait(hCircle2); %confirm position
        innermask = hCircle2.createMask;
        % create annulus outer-inner
        outermask(innermask) = 0;
end

