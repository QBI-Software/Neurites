function [ dcache,FR ] = findCSVdata( T1,Scale,Shiftx,Shifty,pary,xdata,ydata,color )
%Determine csv data within ROI xdata,ydata

        hold on;
%         XC1 = [];
%         YC1 = [];
        FR = [];
        dcache = [];
%         %Border points
%         asize = 1000;
%         px = linspace(x,x+rwidth,asize);
%         py = linspace(y,y+rheight,asize);
%         pleft = [linspace(x,x,asize)',py']; 
%         pright = [linspace(x+rwidth,x+rwidth,asize)',py'];
%         ptop = [px',linspace(y,y,asize)'];
%         pbottom = [px',linspace(y+rheight,y+rheight,asize)'];
%         pary = [pleft pright ptop pbottom];
%         %Show borders
%         plot(pleft(:,1), pleft(:,2),color)
%         plot(pright(:,1), pright(:,2),color)
%         plot(ptop(:,1), ptop(:,2),color)
%         plot(pbottom(:,1), pbottom(:,2),color)
        for j=1: height(T1)
            [xc,yc] = coords2Img(Scale,Shiftx,Shifty,T1.StartX(j),T1.StartY(j));
            [in,on] = inpolygon(xc,yc,xdata,ydata);
            if (in || on)
                plot(xc,yc,'Marker','x','MarkerFaceColor', color)
                FR(end+1) = j;
                if(strcmp(T1.PointType(j),'EP') > 0 )
                   [exc,eyc]= coords2Img(Scale,Shiftx,Shifty,T1.EndX(j),T1.EndY(j));
                    plot(exc,eyc,'Marker','d','MarkerFaceColor', color); %end of neurite  
                else if (~on && j>1)
                    % CALCULATE EXTRA distance at either end as not end 
                    % if previous or next start points are NOT inpolygon
                    % store these minus border
                    extras=[j-1,j+1];
                    for t=1:length(extras)
                        ex = extras(t);
                        [prev_xc, prev_yc]= coords2Img(Scale,Shiftx,Shifty,T1.StartX(ex),T1.StartY(ex)); 
                        if ~inpolygon(prev_xc,prev_yc,xdata,ydata)
                            u = 1;
                            %comp1 = pary(u);
                            pd = pdist2([xc,yc], [pary(:,u), pary(:,u+1)]);
                            [r,c] = find(pd < 100); %10um equiv
                            %u = u+2;
                            [~,maxp] = size(pary); 
                            while (isempty(c) && u < maxp-1)
                                u = u + 2;
                                pd = pdist2([xc,yc], [pary(:,u), pary(:,u+1)]);
                                [r,c] = find(pd < 100);     
                            end
                            if (~isempty(c))
                                if (length(c)>1)
                                    [val,idx] = min(pd(c));
                                    pidx = c(idx);
                                else
                                    pidx = c; %?check
                                end
                                
                                switch t
                                    case 1
                                        dcache(j).begin = [pary(pidx,u),pary(pidx,u+1)];
                                        plot(pary(pidx,u),pary(pidx,u+1),'Marker','+','MarkerFaceColor', color)
                                    case 2
                                        dcache(j).end=[pary(pidx,u),pary(pidx,u+1)];
                                        plot(pary(pidx,u),pary(pidx,u+1),'Marker','o','MarkerFaceColor',color)
                                end
                            end
                        end
                    end
                end
                end %not endpoint
            end %not in or on border

        end

end

