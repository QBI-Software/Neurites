function [ dcache,FR ] = findCSVdata( T1,Scale,Shiftx,Shifty,pary,xdata,ydata,color )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Add csv data
%         T1 = CSVFile;
        
        hold on;
        XC1 = [];
        YC1 = [];
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
        for (j=1: height(T1))
            xc = (T1.StartX(j) * Scale) + Shiftx;
            yc = (T1.StartY(j) * Scale) + Shifty;
            [in,on] = inpolygon(xc,-yc,xdata,ydata);
            if (in || on)
                plot(xc,-yc,'xr')
                XC1(end+1) = xc;
                YC1(end+1) = yc;
                FR(end+1) = j;
                if(strcmp(T1.PointType(j),'EP') > 0 )
                  XC1(end+1) = (T1.EndX(j) * Scale) + Shiftx;
                  YC1(end+1) = (T1.EndY(j) * Scale) + Shiftx;  
                  plot(XC1,-YC1,'color',color,'LineStyle','-','LineWidth', 2);  
                  XC1 = [];
                  YC1 = [];  
                else if (~on && j>1)
                    % CALCULATE EXTRA distance
                    % if previous or next start points are NOT inpolygon
                    % store these minus border
                    prev_xc = (T1.StartX(j-1) * Scale) + Shiftx;
                    prev_yc = (T1.StartY(j-1) * Scale) + Shifty;
                    if ~inpolygon(prev_xc,-prev_yc,xdata,ydata)
                        u = 1;
                        %comp1 = pary(u);
                        pd = pdist2([xc,-yc], [pary(:,1), pary(:,2)]);
                        [r,c] = find(pd < 100); %10um equiv
                        %u = u+2;
                        [~,maxp] = size(pary); 
                        while (isempty(c) && u < maxp-1)
                            u = u + 2;
                            pd = pdist2([xc,-yc], [pary(:,u), pary(:,u+1)]);
                            [r,c] = find(pd < 100);     
                        end
                        if (~isempty(c))
                            if (length(c)>1)
                                [val,idx] = min(pd(c))
                                pidx = c(idx)
                            else
                                pidx = c %?check
                            end
                            plot(pary(pidx,u),pary(pidx,u+1),'xy')
                            cx_x = (pary(pidx,u) - Shiftx)/Scale;
                            cx_y = (pary(pidx,u+1) - Shifty)/Scale;
                            extrapts = [cx_x,cx_y] %should be point on border
                            dcache(j).begin=extrapts;
                        end
                            
                    end
                    %NEXT POINT
                    post_xc = (T1.StartX(j+1) * Scale) + Shiftx;
                    post_yc = (T1.StartY(j+1) * Scale) + Shifty;
                    if ~inpolygon(post_xc,-post_yc,xdata,ydata)
                        u = 1;
                        %comp1 = pary(u);
                        pd = pdist2([xc,-yc], [pary(:,1), pary(:,2)]);
                        [r,c] = find(pd < 100); %10um equiv
                        %u = u+2;
                        [~,maxp] = size(pary); 
                        while (isempty(c) && u < maxp -1)
                            u = u + 2;
                            pd = pdist2([xc,-yc], [pary(:,u), pary(:,u+1)]);
                            [r,c] = find(pd < 100); 
                            
                        end
                        if (~isempty(c))
                            if (length(c)>1)
                                [val,idx] = min(pd(c))
                                pidx = c(idx)
                            else
                                pidx = c %?check
                            end
                            plot(pary(pidx,u),pary(pidx,u+1),'xy')
                            cx_x = (pary(pidx,u) - Shiftx)/Scale;
                            cx_y = (pary(pidx,u+1) - Shifty)/Scale;
                            extrapts = [cx_x,cx_y] %should be point on border
                            dcache(j).end=extrapts;
                        end
                    end
                            
                    end
                end %not endpoint
            end %not in or on border

        end

end

