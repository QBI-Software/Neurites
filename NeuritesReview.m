%% Neurites Review
% Visualize neurites and synaptic locations
function NeuritesReview(I,N,cell1label,cell2label)
    %create image figure
    f = figure('WindowStyle','normal');
    im = imshow(I); 
    hold on;   
    
    i = 1;
    syn = N.Synapses{i};
    plot(N.soma1.centroid(:,1), N.soma1.centroid(:,2),'color', 'b',...
            'marker','o','linestyle','none','LineWidth', 2);
    plot(N.soma2.centroid(:,1), N.soma2.centroid(:,2),'color','b',...
            'marker','o','linestyle','none','LineWidth', 2);
    % Add buttons to gui
    str1 = sprintf('%s soma',cell1label);
    str2 = sprintf('%s end',cell1label);
    str3 = sprintf('%s soma',cell2label);
    str4 = sprintf('%s end',cell2label);
    options = {'Select type',str1,str2,str3,str4};
    c1 = sprintf('%s',cell1label);
    c2 = sprintf('%s',cell2label);
    colouroptions = {'Blue (default)','Red','Yellow','Cyan','Magenta','Black'};

    hp = uipanel('Parent',f,'Title','Review controls','FontSize',12,...
             'Tag','panelReview','BackgroundColor','white',...
             'Units','pixels','Position',[5 5 600 100]);
    tth = uicontrol(hp,'Style','text','String','Select Review Type',...
       'Tag','txtReviewStatus',...
       'Units','pixels','Position',[5 50 100 30]);
    pth = uicontrol(hp,'Style','popup','String',options,...
       'TooltipString','Select type',...
       'Units','pixels','Position',[5 5 100 40],...
       'Callback',{@setreviewtype,syn,N});

    tth1 = uicontrol(hp,'Style','pushbutton','String','Previous',...
       'Units','pixels','Position',[120 50 100 25],...
       'Visible','off','Callback',{@showtrack,-1});
    tth2 = uicontrol(hp,'Style','pushbutton','String','Next',...
       'Units','pixels','Position',[220 50 100 25],...
       'Visible','off','Callback',{@showtrack,1}); 
    tth3 = uicontrol(hp,'Style','checkbox','String','Delete Synapse',...
       'Units','pixels','Position',[400 50 120 25],...
       'Tag','btnReviewDelete',...
       'Visible','off','Callback',{@deleteSynapse});
    tth4 = uicontrol(hp,'Style','pushbutton','String','Remove Branch',...
       'Units','pixels','Position',[230 5 100 25],...
       'Visible','off','Callback',{@removeRegion});

    tth5 = uicontrol(hp,'Style','pushbutton','String','SAVE',...
       'Units','pixels','Position',[330 5 80 25], ...
       'Visible','off','Callback',{@acceptChanges,cell1label, cell2label} );

    tth6 = uicontrol(hp,'Style','pushbutton','String','Clear',...
       'Units','pixels','Position',[320 50 60 25],...
       'Visible','off','Callback',@clearplots );
    tth7 = uicontrol(hp,'Style','pushbutton','String','Show All',...
       'Units','pixels','Position',[500 50 80 25],...
       'Visible','off','Callback',{@showAll} );
      %Select branch for display

    th0 = uicontrol(hp, 'Style','pushbutton','String',strcat('Filter ',c1),...
       'Tag','btnFilter',...
       'Units','pixels','Position',[420 5 80 25],...
       'Visible','off','Callback',{@filterTree,1,c1} );
    th1 = uicontrol(hp, 'Style','pushbutton','String',strcat('Filter ',c2),...
       'Tag','btnFilter',...
       'Units','pixels','Position',[500 5 80 25],...
       'Visible','off','Callback',{@filterTree,2,c2} );
    %Change Soma centroids
    tth8 = uicontrol(hp, 'Style','pushbutton','String','Soma Centroids',...
       'Tag','btnSomaCentroids',...
       'Units','pixels','Position',[120 5 100 25],...
       'Visible','on','Callback',{@changeSoma} );
    %Change colour of overlay
    tth9 = uicontrol(hp,'Style','popup','String',colouroptions,...
       'TooltipString','Select colour',...
       'Units','pixels','Position',[5 -20 100 40],...
       'Callback',{@setcolour});

  
end

function setcolour(source,callbackdata)
    c = source.Value;
    colours = ['b','r','y','c','m','k'];
    hId = findobj('Tag','btnReview');
    hData = get(hId,'UserData');
    hData.colour =colours(c);
    set(hId,'UserData',hData);
end    
function setreviewtype(source,callbackdata,syn,N1)
    c = source.Value -1;
    hId = findobj('Tag','btnReview');
    hData = get(hId,'UserData');
    if (~isempty(hData))
        p = hData.p;
        s = hData.s;
        delete(p);
        delete(s);
        deleted = hData.deleted;
        changed = hData.changed;
        colour = hData.colour;
    else
        deleted = [];
        changed = [];
        colour = 'b';
    end
    linestyle=':';
    [s1,p] = plottrace(c,syn,linestyle,colour);
    if c <=2
        csvfile = N1.CSV1;
    else
        csvfile = N1.CSV2;
    end
    %Review Status
    status = sprintf('Synapse %d of %d (%0.02f um)', 1, length(N1.Synapses),getSynDistance(c,syn));
    hR = findobj('Tag','txtReviewStatus');
    hR.String = status;
    %Enable buttons
    hBtn = findobj('Tag','panelReview');
    kids = allchild(hBtn);
    for i=1:length(kids)
    	kids(i).Visible = 'on';
    end 
    %save data
    reviewdata = struct('i',1,'s',s1,'reviewtype',c,'p',p,'csvfile',csvfile,...
        'deleted',deleted,'changed',changed,'colour',colour);
    set(hId,'UserData',reviewdata);
end
function d = getSynDistance(c,syn)
    switch c
        case 1
            d = syn.SomaC1;
        case 2
            d = syn.NeuriteEndC1;
        case 3
            d = syn.SomaC2;
        case 4
            d = syn.NeuriteEndC2;
    end
end  
function [s1,p] = plottrace(c,syn,linestyle,colour)
    if (isempty(colour))
        colour = 'b';
    end
    if (isempty(linestyle))
        linestyle = ':';
    end
    switch c
        case 1
            s1 = plot(syn.MedianC1(:,1), syn.MedianC1(:,2),'color',...
                'm','marker','X','linestyle','none','LineWidth', 2);
            %synapse may be on soma
            if (isempty(syn.SomapointsC1))
                p = plot(syn.MedianC1(:,1), syn.MedianC1(:,2),...
                'color','y','marker','O','linestyle','none','LineWidth', 1);
            else
                p = plot(syn.SomapointsC1(:,1),syn.SomapointsC1(:,2),...
                    'color',colour,'linestyle',linestyle,'LineWidth', 2);
            end
        case 2
            p = plot(syn.EndpointsC1(:,1),syn.EndpointsC1(:,2),...
                'color',colour,'linestyle',linestyle,'LineWidth', 2);
            s1 = plot(syn.MedianC1(:,1), syn.MedianC1(:,2),'color',...
                'm','marker','X','linestyle','none','LineWidth', 2);
        case 3
            s1 = plot(syn.MedianC2(:,1), syn.MedianC2(:,2),'color',...
                'm','marker','X','linestyle','none','LineWidth', 2);
            %synapse may be on soma
            if (isempty(syn.SomapointsC2))
                p = plot(syn.MedianC2(:,1), syn.MedianC2(:,2),...
                'color','y','marker','O','linestyle','none','LineWidth', 1);
            else
                p = plot(syn.SomapointsC2(:,1),syn.SomapointsC2(:,2), ...
                'color',colour,'linestyle',linestyle,'LineWidth', 2);
            end
            
        case 4
            s1 = plot(syn.MedianC2(:,1), syn.MedianC2(:,2),'color',...
                'm','marker','X','linestyle','none','LineWidth', 2);
            p = plot(syn.EndpointsC2(:,1),syn.EndpointsC2(:,2), ...
                'color',colour,'linestyle',linestyle,'LineWidth', 2);
    end
end     
function p = showtrack(source,callbackdata,ix)
    hId = findobj('Tag','btnIdentify');
    hNData = get(hId,'UserData');
    N = hNData.N;
    hId = findobj('Tag','btnReview');
    hData = get(hId,'UserData');
    i = hData.i;
    if (i >= 1)
        s = hData.s;
        p = hData.p;
        delete(s);
        delete(p);
    end
    i = i + ix;
    %no change if at limits
    if (i < 1 || i > length(N.Synapses))
        i = i - ix; 
    end
    c = hData.reviewtype;
    syn = N.Synapses{i};
    linestyle=':';
    colour = hData.colour;
    [s,p] = plottrace(c,syn,linestyle,colour);
    
    %Review Status
    status = sprintf('Synapse %d of %d (%0.02f um)', i, length(N.Synapses),getSynDistance(c,syn));
    hR = findobj('Tag','txtReviewStatus');
    hR.String = status;
    hD = findobj('Tag','btnReviewDelete');
    if (ismember(i,hData.deleted))
        hD.Value = 1;
        set(p, 'color', [0.5 0.5 0.5]);
    else
        hD.Value = 0;
    end
    hData.p = p;
    hData.i = i;
    hData.s = s;
    
    set(hId,'UserData',hData);
end    
function changeSoma(source,callbackdata)
    hId = findobj('Tag','btnIdentify');
    hNData = get(hId,'UserData');
    N1 = hNData.N;
    %N1 = nH.UserData.analyser;
    s1 = [N1.soma1.centroid(:,1), N1.soma1.centroid(:,2)]
    s2 = [N1.soma2.centroid(:,1), N1.soma2.centroid(:,2)]
    prompt = {'Enter Soma 1 x,y coords:','Enter Soma 2 x,y coords:'};
    dlg_title = 'Change Soma centroid positions';
    num_lines = 1;
    defaultans = {num2str(s1),num2str(s2)};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans)
    %save results
    a1 = strsplit(answer{1},' ');
    a2 = strsplit(answer{2},' ');
    N1.soma1.centroid = [str2num(a1{1}) str2num(a1{2})]
    N1.soma2.centroid = [str2num(a2{1}) str2num(a2{2})]
    clearplots();
    plot(N1.soma1.centroid(:,1), N1.soma1.centroid(:,2),'color', 'b',...
            'marker','o','linestyle','none','LineWidth', 2);
    plot(N1.soma2.centroid(:,1), N1.soma2.centroid(:,2),'color','b',...
            'marker','o','linestyle','none','LineWidth', 2);
    %save back
    N1 = N1.updateSomaCalculations();
    hNData.N = N1;
    set(hId,'UserData',hNData);
end    


function filterTree(source,callbackdata,cellnum, celllabel)
    hNId = findobj('Tag','btnIdentify');
    hNData = get(hNId,'UserData');
    N1 = hNData.N;
    branchoptions = { N1.getbranchoptions(1) N1.getbranchoptions(2)};   
    prompt = sprintf('Select %s tree-branch:',celllabel);
    [s,v] = listdlg('PromptString',prompt,...
                'SelectionMode','multiple',...
                'CancelString', 'None',...
                'ListString',branchoptions{cellnum})
    %save filters
    if (v) 
        b = branchoptions{cellnum}(s);
        N1 = N1.setfilter(cellnum,b);
    
        %show selected tree-branches
        hId = findobj('Tag','btnReview');
        hData = get(hId,'UserData');
        colour = hData.colour;
        clearplots();
        for i=1:length(N1.Synapses)
            syn = N1.Synapses{i};
            if cellnum == 1
                tb = strcat(num2str(syn.TreeC1),'-',num2str(syn.BranchPointC1));
                c = 1
            else
                tb = strcat(num2str(syn.TreeC2),'-',num2str(syn.BranchPointC2));
                c = 3;
            end
            if (ismember(tb,b))
                [s,p] = plottrace(c,syn,'-',colour);
            end
        end
    else
        %clear filters
        b = {}
        N1 = N1.setfilter(cellnum,b);
    end
    %save
    hNData.N = N1;
    set(hNId,'UserData',hNData);
end    
function showAll(source,callbackdata)
    hId = findobj('Tag','btnReview');
    hData = get(hId,'UserData');
    hId = findobj('Tag','btnIdentify');
    hNData = get(hId,'UserData');
    N1 = hNData.N;
    c = hData.reviewtype;
    colour = hData.colour;
    clearplots();
    plot(N1.soma1.centroid(:,1), N1.soma1.centroid(:,2),'color', 'b',...
            'marker','o','linestyle','none','LineWidth', 2);
    plot(N1.soma2.centroid(:,1), N1.soma2.centroid(:,2),'color','b',...
            'marker','o','linestyle','none','LineWidth', 2);
    for i=1:length(N1.Synapses)
        syn = N1.Synapses{i};
        [s,p] = plottrace(c,syn,'-',colour);
         if (ismember(i,hData.deleted))
            set(p, 'color', [0.5 0.5 0.5]);
         end
    end
    %Also show ROI
    hRoi = findobj('Tag','radioROI');
    %roidata = hRoi.UserData;
    roidata = get(hRoi,'UserData');
    if (~isempty(roidata))
        xi = roidata.xi;
        yi = roidata.yi;
        plot(xi, yi, '--b','LineWidth', 1);
    end
end    
function clearplots(source,callbackdata)
    cla;
    h = findobj('Tag','btnBrowser');
    hData = get(h,'UserData');
    I = hData.img;
    %f = figure('WindowStyle','normal');
    im = imshow(I);
end 
%Mark for deletion (allows reversible)
function deleteSynapse(source, callbackdata)
    hId = findobj('Tag','btnReview');
    hData = get(hId,'UserData');
    i = hData.i;
    p = hData.p;
    s = hData.s;
    %check toggle value
    if (source.Value)
        hData.deleted(end+1) = i;
        set(p, 'color', [0.5 0.5 0.5]);
        set(s, 'color', [0.5 0.5 0.5]);
    else
        remove = find(ismember(hData.deleted,i));
        if (remove > 0)
            hData.deleted(remove) = [];
            set(p, 'color', 'b');
            set(s, 'color', 'm');
        end
    end
            
    set(hId,'UserData',hData);
end   
function N1 = removeRegion(source,callbackdata)
    hId = findobj('Tag','btnIdentify');
    hNData = get(hId,'UserData');
    N1 = hNData.N;
    hRId = findobj('Tag','btnReview');
    hData = get(hRId,'UserData');
    i = hData.i;
    c = hData.reviewtype;
    p = hData.p;
    s = hData.s;
    syn = N1.Synapses{i};
    csvfile = hData.csvfile;
    dcm_obj = datacursormode(gcf);
    set(dcm_obj,'DisplayStyle','datatip',...
        'SnapToDataVertex','on','Enable','on')
    msgbox('Click on a branch to remove it, then click Enter.')
    % Wait while the user does this.
    pause 
    c_info = getCursorInfo(dcm_obj);
    cachesyn = syn;
    syn = syn.removeBranch(c,csvfile,...
      c_info.Target.XData(c_info.DataIndex),...
      c_info.Target.YData(c_info.DataIndex));
    delete(p);
    delete(s);
    linestyle=':';
    [s,p] = plottrace(c,syn,linestyle);
    
    v = accept(i);
    if(v==1)
        dcm_obj.Enable='off';
        N1.Synapses{i} = syn; %update
        hData.changed(end+1) = i;
        hNData.N = N1;
        set(hId,'UserData',hNData);
    else %reject change
         syn = cachesyn;
         [s,p] = plottrace(c,syn,linestyle);
    end
    hData.s = s;
    hData.p = p;
    set(hRId,'UserData',hData);
end
function a = accept(synnum)
    prompt = sprintf('Accept this measurement for synapse %d?',synnum);
    %str = input(prompt,'s');
    str = questdlg(prompt,'Synapse change',...
        'Yes','No','Yes');
    switch str
        case 'Yes'
            a = 1;
        case 'No'
            a = 0;
    end
end
function acceptChanges(source,callbackdata,cell1label, cell2label)
    prompt = 'Accept changes and update results (with filters)?';
    str = questdlg(prompt,'Update results',...
        'Yes','No','Yes');
    hId = findobj('Tag','btnReview');
    hData = get(hId,'UserData');
    
    switch str
        case 'Yes'
             hNId = findobj('Tag','btnIdentify');
             hNData = get(hNId,'UserData');
             N1 = hNData.N;
             %remove deleted
             for d=1:length(hData.deleted)
                 i = hData.deleted(d);
                 N1.Synapses{i} = [];
             end
             [colnames,T] = N1.generateTable([1],cell1label, cell2label);    
             htable = findobj('Tag','uitableResults');
             set(htable,'data',T,'ColumnName',colnames);
             %save to file
             hCSV = findobj('Tag','btnAnalysisFiles');
             csvdata = get(hCSV,'UserData');
             pathname = csvdata.csvPath;
             outputfile = fullfile(pathname, 'neurites_data_review.csv');
             [FileName,PathName] = uiputfile({'*.csv','CSV file'},...
                 'Save Data', outputfile)
             outputfile = fullfile(PathName,FileName);
             saveDataFile(outputfile, colnames, T);
             %Recopy synapses to new list
             if (length(hData.deleted))
                 cSynapses ={length(N1.Synapses)};
                 m = 1;
                 for j=1:length(N1.Synapses)
                     syn = N1.Synapses{j};
                     if (~isempty(syn))
                         cSynapses{m} = syn;
                         m = m+1;
                     end
                 end
                 hNData.N.Synapses = cSynapses;
             else
                 hNData.N = N1;
             end
             
             msg = sprintf('Changed %d synaptic regions. Deleted %d synaptic regions. \nData updated to %s', length(hData.changed),length(hData.deleted), outputfile);
             msgbox(msg);
             set(hNId,'UserData',hNData);
             %Reset deleted and changed
             hData.i = 1;
             hData.deleted=[];
             hData.changed=[];
             set(hId,'UserData',hData);
             
             close(source.Parent.Parent); %close figure
        case 'No'
            a = 0;
    end
end