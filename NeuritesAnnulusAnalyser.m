classdef NeuritesAnnulusAnalyser
   properties
      I
      Iroi
      bw
      maskedI
      CSV
      soma
      Regions
      filter1
      centreshiftx
      centreshifty
   end
   methods
      %Example: N = NeuritesAnalyser(image_as_rgb);
      function obj = NeuritesAnnulusAnalyser(rgbimg)
         if (ischar(rgbimg))
          obj.I = imread(rgbimg);
         else
          obj.I = rgbimg;
         end
         
         % create binary version of image
         obj.bw = im2bw(obj.I);
          
         % Find cell somas
         obj.soma = obj.findSoma(obj.bw);
         
         %Unload from memory
         obj.I = 0;
      end
      
      % ROI mask can be read from a file or passed directly
      % If annulus - create mask from inner and outer circles first
      function obj = loadMask(obj, roi)
          
         if (ischar(roi)) 
          obj.Iroi = imread(roi);
         else
          obj.Iroi = roi;
         end
         % Create masked image from roi
         obj.maskedI = obj.bw; 
         obj.maskedI(~obj.Iroi) = 0;  % Set all non-keeper pixels to zero.
      end
      
      
      
      function soma = findSoma(obj, rbw)
        %Isolate soma for analysis by removing dendrites
        %BW = imfill(rbw,'holes');
        BW = bwmorph(rbw,'majority');
        se1 = strel('rectangle',[10 2]);
        se2 = strel('rectangle',[2 10]);
        se3 = strel('square',30);
        BW1 = imopen(BW,se1);
        BW1 = imopen(BW1,se2);
        [B,L,N,~] = bwboundaries(BW1);%,8,'noholes');
        if (N > 30) %thick dendrites - be brutal
            BW1 = imopen(BW1,se3);
            [B,L,N,~] = bwboundaries(BW1);
        end
        
        stats=  regionprops(L, 'Centroid', 'Area', 'Perimeter',...
            'MajorAxisLength','MinorAxisLength','Eccentricity','Solidity');
        Centroid = cat(1,stats.Centroid);
        Perimeter = cat(1,stats.Perimeter);
        Eccentricity = cat(1,stats.Eccentricity);
        Area = cat(1,stats.Area);
        %Radii = cat(1,sqrt(Area/pi));
        %circles = find(Radii > 5 & Radii < 100);
        idx = find(Area == max(Area));
        %Iterate to find another suitable region  TODO
        B = sort(Area,'descend');
        i = 2;
        while (Eccentricity(idx) > 0.8 && i<= length(Area) ) %not a circle
            %sort Area values by max first then take second
            idx = find(Area == B(i));
            i = i+1;
        end
               
        c = Centroid(idx,:);
        
        soma = NeuritesSoma(Area(idx),Perimeter(idx),c)
%         figure
%         imshow(label2rgb(L, @jet, [.5 .5 .5]))
%         hold on
%         plot(soma.centroid(:,1), soma.centroid(:,2),'color', 'r',...
%             'marker','o','linestyle','none','LineWidth', 2);
        
      end
      
      function obj = set.Regions(obj,segments)
          if (~isempty(segments))
             obj.Regions = segments;
          else
             error('Error: Cannot set REGIONS - Empty set')
          end
      end
    
    %% Generate Table data - export data to a table format
    function [colnames,tabledata] = generateTable(obj,cell1label, cell2label)
         colnames = {'Cell1_X' 'Cell1_Y' 'Cell2_X' 'Cell2_Y',...
             'Cell1_tree' 'Cell1_order' 'Cell1_length' 'Cell1_distance' 'Cell1_soma' 'Cell1_end'  ...
             'Cell1_theta' 'Cell1_rho' 'Cell1_deg' ...
             'Cell2_tree' 'Cell2_order' 'Cell2_length' 'Cell2_distance' 'Cell2_soma' 'Cell2_end'  ...
             'Cell2_theta' 'Cell2_rho' 'Cell2_deg' ...
             'Cell1_StartX' 'Cell1_StartY' 'Cell2_StartX' 'Cell2_StartY' };
         colnames = strrep(colnames, 'Cell1', cell1label);
         colnames = strrep(colnames, 'Cell2', cell2label);
%         colnames = {'Type' 'DS_X' 'DS_Y' 'SBAC_X' 'SBAC_Y' ...
%             'DS_length' 'DS_distance' 'DS_order' 'DS_end' 'DS_soma'  ...
%             'SBAC_length' 'SBAC_distance' 'SBAC_order' 'SBAC_end' 'SBAC_soma' ...
%             'DS_StartX' 'DS_StartY' 'SBAC_StartX' 'SBAC_StartY' };
        tabledata =zeros(length(obj.Regions),length(colnames)); %preallocate
        ignorefilter=0;
        %Override filters if empty
        if (isempty(obj.filter1) && isempty(obj.filter2))
            ignorefilter=1;
        end
        for i=1:length(obj.Regions)
            syn = obj.Regions{i};
            if (~isempty(syn))
                %Apply filters
                tb1 = strcat(num2str(syn.TreeC1),'-',num2str(syn.BranchPointC1));
                tb2 = strcat(num2str(syn.TreeC2),'-',num2str(syn.BranchPointC2));
                if (ignorefilter || (~isempty(obj.filter1) && ismember(tb1,obj.filter1)) ...
                        || (~isempty(obj.filter2) && ismember(tb2,obj.filter2)))
                    %show data in table
                    [syn1coordsx, syn1coordsy]  = syn.img2Coords(syn.MedianC1(1),syn.MedianC1(2));
                    [syn2coordsx, syn2coordsy] = syn.img2Coords(syn.MedianC2(1),syn.MedianC2(2));
                    row1 = [syn1coordsx syn1coordsy syn2coordsx syn2coordsy ...
                        syn.TreeC1 syn.BranchPointC1 syn.NeuriteLengthC1 syn.DistanceC1 ...
                        syn.SomaC1 syn.NeuriteEndC1 ...
                        syn.ThetaC1 syn.RhoC1 syn.DegC1 ...
                        syn.TreeC2 syn.BranchPointC2 syn.NeuriteLengthC2 syn.DistanceC2 ...
                        syn.SomaC2 syn.NeuriteEndC2 ...
                        syn.ThetaC2 syn.RhoC2 syn.DegC2 ...
                        syn.StartXY1 syn.StartXY2];
                    %tabledata = cat(1,tabledata,row1);
                    tabledata(i,:) = row1;
                end
            end
        end
        tabledata( ~any(tabledata,2), : ) = [];  %remove zero rows
    end
    
   
   
   end
end







      
