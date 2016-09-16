function [ ep ] = findEndpoints( bw_image )
%Determine end points of a region via extrema
s = regionprops(bw_image,'Extrema');
plot(s.Extrema(:,1),s.Extrema(:,2),'xr','LineWidth',2)
eps = [];
for i=1:length(s.Extrema)
    a = [s.Extrema(1,i),s.Extrema(2,i)];
    %check next in list
    [D,i] = knnsearch(s.Extrema,a,'K',2,'Distance','euclidean');
    %skip next if too close
    if (i(2) < 20)
        i = i+1;
    end
    eps(i) = a;
end

