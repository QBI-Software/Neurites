function varargout = NeuritesAnnulusAppUI(varargin)
% NEURITESANNULUSAPPUI MATLAB code for NeuritesAnnulusAppUI.fig
%      NEURITESANNULUSAPPUI, by itself, creates a new NEURITESANNULUSAPPUI or raises the existing
%      singleton*.
%
%      H = NEURITESANNULUSAPPUI returns the handle to a new NEURITESANNULUSAPPUI or the handle to
%      the existing singleton*.
%
%      NEURITESANNULUSAPPUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NEURITESANNULUSAPPUI.M with the given input arguments.
%
%      NEURITESANNULUSAPPUI('Property','Value',...) creates a new NEURITESANNULUSAPPUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NeuritesAnnulusAppUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NeuritesAnnulusAppUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% ***NB***% For R2014a and earlier: 
%       data = get(h,'UserData');
%       set(hObject,'UserData',data); 
% Edit the above text to modify the response to help NeuritesAnnulusAppUI

% Last Modified by GUIDE v2.5 09-Jun-2016 14:04:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NeuritesAnnulusAppUI_OpeningFcn, ...
                   'gui_OutputFcn',  @NeuritesAnnulusAppUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before NeuritesAnnulusAppUI is made visible.
function NeuritesAnnulusAppUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NeuritesAnnulusAppUI (see VARARGIN)

% Choose default command line output for NeuritesAnnulusAppUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using NeuritesAnnulusAppUI.
if strcmp(get(hObject,'Visible'),'off')
    %(rand(5));
    %load image in panel
    axes(handles.axes1);
    hFig = figure(get(handles.axes1,'parent'));%handles.axes1.Parent;
    hIm = imshow('synapse_EM.jpg');
    hSP = imscrollpanel(hFig,hIm); % Handle to scroll panel.
    set(hSP,'Units','pixels',...
        'Position',[50 190 551 450]) %match with GUI box
	    
    % Add a Magnification Box and an Overview tool.
    hMagBox = immagbox(hFig,hIm);
    pos = get(hMagBox,'Position');
    set(hMagBox,'Position',[0 0 pos(3) pos(4)])
    imoverview(hIm)
    %Get the scroll panel API to programmatically control the view.
    api = iptgetapi(hSP);
    %Get the current magnification and position.
    %mag = api.getMagnification();
    %r = api.getVisibleImageRect();
    %View the top left corner of the image.
    %api.setVisibleLocation(0.5,0.5)
    %Change the magnification to the value that just fits.
    api.setMagnification(api.findFitMag())
    %Zoom in to 1600% on the dark spot.
    %api.setMagnificationAndCenter(16,306,800)
    %test
    I = imread('synapse_EM.jpg');
    api.replaceImage(I)
    
end

% UIWAIT makes NeuritesAnnulusAppUI wait for user response (see UIRESUME)
 %uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = NeuritesAnnulusAppUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
%function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%file = uigetfile('*.fig');
%if ~isequal(file, 0)
%    open(file);
%end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

%set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});

% --------------------------------------------------------------------
function menu_File_loadimage_Callback(hObject, eventdata, handles)
% hObject    handle to menu_File_loadimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    fileTypes = {  '*.tif;*.tiff;*.jpg;*.jpeg', 'Tiff/JPG images' };
    
    [imageFile, imagePath, ~] = ...
      uigetfile(fileTypes, ...
		'Select neuron image', ...
		'MultiSelect', 'off');
    inputPath = fullfile(imagePath,imageFile);
    updateStatus(handles,'Loading image ....');
    I = imread(inputPath);
    %load image in panel
    axes(handles.axes1);
    hSP = get(handles.axes1,'parent');
    api = iptgetapi(hSP);
    api.replaceImage(I);
    % Reset magnification
    mag = api.findFitMag()
    api.setMagnification(mag);
    htable = findobj('Tag','uitableResults');
    set(htable,'Data',[]);
    %Clear configuration
    clearConfig(handles,0);
    clearCSVdata(handles);
    %Directions
    updateStatus(handles,'Image loaded.');
    %Save image analyser
    N = NeuritesAnnulusAnalyser(I); 
    %save image data
    data = struct('img',I, 'imagePath',imagePath,'inputfile', inputPath, 'analyser', N);
    set(hObject,'UserData',data); 

% --------------------------------------------------------------------
function Menu_File_loadcsv_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_File_loadcsv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    fileTypes = {  '*.csv', 'CSV' };
    hfiles = findobj('Tag','menu_File_loadimage');
    %hfdata = hfiles.UserData;
    hfdata = get(hfiles,'UserData');
    csvPath = hfdata.imagePath;
    [csvFile, csvPath, ~] = ...
      uigetfile(fileTypes, ...
		'Select analysis file', ...
		'MultiSelect', 'off',csvPath);
    csv1 = fullfile(csvPath, csvFile);
    
    %Check CSV headers
    hdrs = 'Tree,Order,StartX,StartY,StartZ,EndX,EndY,EndZ,PointType,Length(µm),LengthToBeginning(µm)';
    r1 = checkHeaders(csv1,hdrs);
    if (length(r1) > 1)
        errhdr = sprintf('ERROR: CSV file headers need to match: %s. \nFields incorrect:%s %s', hdrs,r1);
        errordlg(errhdr,'CSV Files error')
        return 
    end
    %save to userdata
    data = get(hObject,'UserData');
    if (isempty(data))
        data = struct('csvPath', csvPath, 'csvFile', csvFile);
    else
        data.csvPath = csvPath;
        data.csvFile = csvFile;
   
    end
    %hObject.UserData = data;
    set(hObject,'UserData',data);
    if (isfield(data,'csvFile'))
         msgbox('CSV files loaded. Register image with CSV data with Register CSV.')
         status = sprintf('CSV loaded: %s', csv1);
         updateStatus(handles,status);
    else
        msgbox('Error in loading CSV file')
    end


% --- Executes on button press in btnRegister.
function btnRegister_Callback(hObject, eventdata, handles)
% hObject    handle to btnRegister (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = findobj('Tag','menu_File_loadimage');
%data = h.UserData;
data = get(h,'UserData');
I = data.img;
N = data.analyser;
%csvdata = handles.btnAnalysisFiles.UserData;
hCSV = findobj('Tag','Menu_File_loadcsv');
csvdata = get(hCSV,'UserData');

if ( ~isempty(csvdata))
    csv1 = fullfile(csvdata.csvPath, csvdata.csvFile);
    %Read values
    hScale = findobj('Tag','editScale');
    hSX = findobj('Tag','editShiftx');
    hSY = findobj('Tag','editShifty');
    
    if isempty(get(hScale,'String'))
        %load config data if empty
        configfile =fullfile(csvdata.csvPath, 'neurites_annulus_config.csv');
        if (exist(configfile, 'file') == 2)
            M = readtable(configfile);
            N.soma.centroid = [M.CentroidX M.CentroidY]; %update centroid
            data.analyser = N;
            set(h,'UserData',data);
        else
            M = estimateOverlay(handles,I,csv1);
            if (exists(data.analyser))
                M.CentroidX = N.soma.centroid(1);
                M.CentroidY = N.soma.centroid(2);
            end
        end
        loadConfig(M,handles);
        
        hscale = M.Scale;
        shiftx = M.Shiftx;
        shifty = M.Shifty;
       
    else
        hscale = str2double(get(hScale,'String'));
        shiftx = str2double(get(hSX,'String'));
        shifty = str2double(get(hSY,'String'));
        
    end  
    
    T1 = readtable(csv1);
    %plot
    axes(handles.axes1); 
    hIm = imshow(I);
    %imoverview(hIm)
    hold on;
    XC1 = [];
    YC1 = [];
    for (i=1: height(T1))
        XC1(end+1) = (T1.StartX(i) * hscale) + shiftx;
        YC1(end+1) = (T1.StartY(i) * hscale) + shifty;
        
        if(strcmp(T1.PointType(i),'EP') > 0)
          plot(XC1,-YC1,'color','c','LineStyle','-','LineWidth', 2);  
          XC1 = [];
          YC1 = [];  
        end
        
    end
    hold off;
end
   
 
    
%Estimate parameters for overlay of CSV data onto image
function M = estimateOverlay(handles, I,cell1csv)
    Img = rgb2gray(I);
    Tb = readtable(cell1csv);
    
    %CSV data in um
    mi = min(Tb.StartX);
    if (min(Tb.EndX) < mi)
        rowi = find(Tb.EndX == min(Tb.EndX));
        S1 = [Tb.EndX(rowi),Tb.EndY(rowi)];
    else
        rowi = find(Tb.StartX == min(Tb.StartX));
        S1 = [Tb.StartX(rowi),Tb.StartY(rowi)];
    end
    mx = max(Tb.StartX);
    if (max(Tb.EndX) > mx)
        rowi = find(Tb.EndX == max(Tb.EndX));
        S2 = [Tb.EndX(rowi),Tb.EndY(rowi)];
    else
        rowi = find(Tb.StartX == max(Tb.StartX));
        S2 = [Tb.StartX(rowi),Tb.StartY(rowi)];
        
    end
    if (length(S1(:,1)) > 1)
        rowi = find(S1(:,2) == max(S1(:,2)));
        S1 = [S1(rowi,1),S1(rowi,2)];
    end
    if (length(S2(:,1)) > 1)
        rowi = find(S2(:,2) == max(S2(:,2)));
        S2 = [S2(rowi,1),S2(rowi,2)];
    end
    C1 = cat(1,S1,S2);
    csvlength = pdist(C1,'euclidean');
    %Image coords - 
    LI = bwlabel(Img);
    s = regionprops(LI, 'Area', 'BoundingBox','Extrema');
    i = 1;
    Area = cat(1,s.Area);
    if (length(s) > 1)
        i = find(Area==max(Area));
    end
%     S3=s(i).Extrema(1,:);
%     S4=s(i).Extrema(4,:);
%     S3=[s(i).BoundingBox(1) s(i).BoundingBox(2)];
%     S4=[s(i).BoundingBox(1)+s(i).BoundingBox(3)  s(i).BoundingBox(2)+s(i).BoundingBox(4)];
    Ic = corner(Img,'Harris','SensitivityFactor',0.02);
    xi = Ic(:,1)==min(Ic(:,1));
    S3 = Ic(xi,:)
    xm = Ic(:,1)==max(Ic(:,1));
    S4 = Ic(xm,:)
    if (length(S3(:,1)) > 1)
        rowi = find(S3(:,2) == max(S3(:,2)));
        S3 = [S3(rowi,1),S3(rowi,2)];
    end
    if (length(S4(:,1)) > 1)
        rowi = find(S4(:,2) == max(S4(:,2)));
        S4 = [S4(rowi,1),S4(rowi,2)];
    end
    % Find scale
    C2 = cat(1,S3,S4);
    imglength = pdist(C2,'euclidean');
    Scale = round(imglength/csvlength,1)
    S = S1 * Scale;
    
    %Plot
    axes(handles.axes1);
    imshow(Img);
    hold on;
    %plot(S(1),S(2),'color','b','Marker','o','LineWidth', 2);
    plot(S3(1),S3(2),'color','y','Marker','x','LineWidth', 2);
    plot(S4(1),S4(2),'color','r','Marker','x','LineWidth', 2);
    hold off;
    %Load data
    Fit = 10;
    
    Shiftx = round(S3(1) - S(1),2);
    Shifty = round(-S3(2) - S(2),2);
    
    
    M = table(Scale, Fit, Shiftx, Shifty);
    
        
        
function loadConfig(M, handles)
    hScale = findobj('Tag','editScale');
    hSX = findobj('Tag','editShiftx');
    hSY = findobj('Tag','editShifty');
    hFit = findobj('Tag','editFit');
    hCx = findobj('Tag','editCentroidX');
    hCy = findobj('Tag','editCentroidY');
    hOD = findobj('Tag','editOD');
    hID = findobj('Tag','editID');
    set(hScale,'String',num2str(M.Scale));
    set(hFit,'String',num2str(M.Fit));
    set(hSX,'String',num2str(M.Shiftx));
    set(hSY,'String',num2str(M.Shifty));
    set(hCx,'String',num2str(M.CentroidX));
    set(hCy,'String',num2str(M.CentroidY));
    set(hOD,'String',num2str(M.AnnulusOD));
    set(hID,'String',num2str(M.AnnulusID));
    

function clearConfig(handles,clearlabels)
    initval = '';
    hScale = findobj('Tag','editScale');
    hSX = findobj('Tag','editShiftx');
    hSY = findobj('Tag','editShifty');
    hFit = findobj('Tag','editFit');
    hCx = findobj('Tag','editCentroidX');
    hCy = findobj('Tag','editCentroidY');
    hOD = findobj('Tag','editOD');
    hID = findobj('Tag','editID');
    
    set(hScale,'String',initval);
    set(hFit,'String',initval);
    set(hSX,'String',initval);
    set(hSY,'String',initval);
    set(hCx,'String',initval);
    set(hCy,'String',initval);
    set(hOD,'String',initval);
    set(hID,'String',initval);

function clearCSVdata(handles)
    hCf0 = findobj('Tag', 'Menu_File_loadcsv');
    set(hCf0,'UserData',{});
    

    % --- Executes on button press in pushbutton9.
function saveConfig_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    %Results table
    hScale = findobj('Tag','editScale');
    hSX = findobj('Tag','editShiftx');
    hSY = findobj('Tag','editShifty');
    hFit = findobj('Tag','editFit');
    hCx = findobj('Tag','editCentroidX');
    hCy = findobj('Tag','editCentroidY');
    hOD = findobj('Tag','editOD');
    hID = findobj('Tag','editID');
    Scale = str2double(get(hScale,'String'));
    Fit = str2double(get(hFit,'String'));
    Shiftx = str2double(get(hSX,'String'));
    Shifty = str2double(get(hSY,'String'));
    CentroidX = str2double(get(hCx,'String'));
    CentroidY = str2double(get(hCy,'String'));
    AnnulusOD = str2double(get(hOD,'String'));
    AnnulusID = str2double(get(hID,'String'));

    T = table(Scale, Fit, Shiftx, Shifty, CentroidX, CentroidY, AnnulusOD, AnnulusID);
    %save to file
    %csvdata = handles.btnAnalysisFiles.UserData;
    hCSV = findobj('Tag','Menu_File_loadcsv');
    csvdata = get(hCSV,'UserData');
    outputfile = 'neurites_annulus_config.csv';
    if ( ~isempty(csvdata))
        outputfile =fullfile(csvdata.csvPath, outputfile);
    end
    writetable(T,outputfile);
    status = sprintf('Config file saved:%s',outputfile);
    updateStatus(handles,status);
    
function [cell1,cell2]=generateHistogram(I, handles)
    %Generate histogram
    IG = rgb2gray(I);
    axes(handles.axes2);
    cla;
    
    %histogram - exclude white
    [px,levels] = imhist(IG);
    bar(px);
    %histogram - exclude white - this should be configurable?
    xlim([0 254]);
    xlabel('color values')
    ylabel('pixels')
    grid on;
    
    title('Image Histogram');
    hBg = findobj('Tag','radioWhite');
    bg = 256;
    if (get(hBg,'Value') == 0)
        bg = 0;
    end
    [cell1,cell2] = SplitCells(I,bg);


function rtn = checkHeaders(csvfile,hdrs)
    rtn = 0;           
    fid = fopen(csvfile);
    C = textscan(fid, '%s %s %s %s %s %s %s %s %s %s %s', 'delimiter', ',', ...
                 'treatAsEmpty', {'NA', 'na'}, ...
                 'commentStyle', '//');
    fclose(fid);
    for i=1: length(C)
        v = strjoin(C{1,i}(1));
        v = strrep(v,' ','');
        if (size(strfind(hdrs,v)) == 0)
            csvfile
            rtn = v
            break
        end
    end
    
        

% --------------------------------------------------------------------

% --- Executes on button press in radioZoom.
function radioZoom_Callback(hObject, eventdata, handles)
% hObject    handle to radioZoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioZoom
if (get(hObject,'Value') > 0)
    zoom on
else
    zoom off
end


% --- Executes on button press in btnImtool.
function btnImtool_Callback(hObject, eventdata, handles)
% hObject    handle to btnImtool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = findobj('Tag','menu_File_loadimage');
%data = h.UserData;
data = get(h,'UserData');
imtool(data.img)


% --- Executes on button press in radioPan.
function radioPan_Callback(hObject, eventdata, handles)
% hObject    handle to radioPan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioPan
if (get(hObject,'Value') > 0)
    pan on
else
    pan off
end
% --- Executes on button press in radioDatacursor.
function radioDatacursor_Callback(hObject, eventdata, handles)
% hObject    handle to radioDatacursor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioDatacursor
if (get(hObject,'Value') > 0)
    datacursormode on
else
    datacursormode off
end
%dcm_obj = datacursormode(handles.axes1)

% --- Executes on button press in radioROI.
function radioROI_Callback(hObject, eventdata, handles)
% hObject    handle to radioROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioROI
if (get(hObject,'Value') > 0)
    axes(handles.axes1);
    %set interactive polygon tool
    [x,y, BW, xi, yi] = roipoly;
    
    %ROI stats
    numPixels = sum(BW(:))
    roi_area = polyarea(xi,yi);
    %temp - copy to file
    hfiles = findobj('Tag','btnBrowser');
    %hfdata = hfiles.UserData;
    hfdata = get(hfiles,'UserData');
    roifile = fullfile(hfdata.imagePath, 'neurites_roi.tif');
    imwrite(BW, roifile);
    %save data
    dataroi = struct('roi', BW, 'xi',xi,'yi',yi,'x',x,'y',y);
    %hObject.UserData = dataroi;
    set(hObject,'UserData',dataroi);
    %show mask image in plot 2
    %setTitlePlot2(handles,'');   
    axes(handles.axes2);
    cla;
    imshow(BW)
    title('ROI Mask');
    %Convert ROI to um2 if csv config loaded
    hScale = findobj('Tag', 'editScale');
    scale = str2double(get(hScale,'String'));
    if (isempty(scale) || scale ==0)
        scale = 1;
    end
    roi_um2 = roi_area / scale;
    %Results table
    RCols={'ROI Number Pixels', 'ROI Area (px)','ROI Area (um2)'};
    RData=[numPixels, roi_area,roi_um2];
    htable = findobj('Tag','uitableResults');
    set(htable,'data',RData,'ColumnName',RCols);
    %Directions
    status= sprintf('ROI set: Proceed to Identify button')
    updateStatus(handles,status);
else
    %hObject.UserData = {};
    set(hObject,'UserData',{});
end

% --- Executes on button press in radioNone.
function radioNone_Callback(hObject, eventdata, handles)
% hObject    handle to radioNone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioNone
if (get(hObject,'Value') > 0)
    zoom off
    pan off
else
    %nothing
end

% --- Executes on button press in radioDraw.
function radioDraw_Callback(hObject, eventdata, handles)
% hObject    handle to radioDraw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioDraw

function outermask = applyAnnulus(handles, id, od, scale, centroid)
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
        updateStatus(handles,status);
        disp(status)
        position = wait(hCircle1); %confirm position
        outermask = hCircle1.createMask;
        %[xmin ymin width height] = position
        %[3995.90774907749 4541.9667896679 3500 3500]
        x = centroid(1) - (rid/2)
        y = centroid(2) - (rid/2)
        inner = [x y rid rid]
        hCircle2 = imellipse(gca, inner);
        status = sprintf('Double-Click on INNER annulus to accept position.');
        updateStatus(handles,status);
        disp(status)
        position = wait(hCircle2); %confirm position
        innermask = hCircle2.createMask;
        % create annulus outer-inner
        outermask(innermask) = 0;
        
       

% --- Executes on button press in radioAuto.
function radioAuto_Callback(hObject, eventdata, handles)
% hObject    handle to radioAuto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioAuto
if (get(hObject,'Value') > 0)
    hOD = findobj('Tag','editOD');
    od = str2double(get(hOD,'String'));
    hID = findobj('Tag','editID');
    id = str2double(get(hID,'String'));
    h = findobj('Tag','menu_File_loadimage');
    data = get(h,'UserData');
    hScale = findobj('Tag','editScale');
    scale = str2double(get(hScale,'String'));
     
    if (isnan(od) || isnan(id))
        msgbox('Error: Please enter OD and ID values first');
    else
        
        N = data.analyser;
        centroid =N.soma.centroid;
        mask = applyAnnulus(handles, id, od, scale, centroid);
        N = N.loadMask(mask);
        imshow(N.maskedI);
        roifile = fullfile(data.imagePath, 'neurites_annulus.tif');
        imwrite(mask, roifile);
        status = sprintf('Annulus mask created: %s.', roifile);
        updateStatus(handles,status);
        disp(status)
        
    end
end

% --- Executes on button press in radioAnnNone.
function radioAnnNone_Callback(hObject, eventdata, handles)
% hObject    handle to radioAnnNone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioAnnNone






function updateStatus(handles,string)
    if strfind(string,'ERROR')
        set(handles.textOutput, 'string', string, 'Foreground', [1 0 0])
    else
        set(handles.textOutput, 'string', string, 'Foreground', [0 0 0])
    end


% --- Executes on button press in btnIdentify.
function btnIdentify_Callback(hObject, eventdata, handles)
% hObject    handle to btnIdentify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clearplots();
% Get grayscale img
h = findobj('Tag','btnBrowser');
%data = h.UserData;
data = get(h,'UserData');
hRoi = findobj('Tag','radioROI');
%roidata = hRoi.UserData;
roidata = get(hRoi,'UserData');
%csvdata = handles.btnAnalysisFiles.UserData;
hCSV = findobj('Tag','btnAnalysisFiles');
csvdata = get(hCSV,'UserData');
if (~isempty(roidata))
    roi = roidata.roi;
    xi = roidata.xi;
    yi = roidata.yi;
else
    if(~isempty(csvdata))
        roi = fullfile(csvdata.csvPath, 'neurites_roi.tif');
        if(~exist(roi, 'file'))
            errordlg('Please set ROI first','ROI error!')
            return
        end
    else
        errordlg('Please set ROI first','ROI error!')
        return
    end
end

% Get configuration
hScale = findobj('Tag','editScale');
if (isempty(get(hScale,'String')))
    errordlg('Please set Configuration first with Register CSV','Configuration error!')
    return
end
h1 = findobj('Tag', 'editMinlength');
minlength = str2double(get(h1,'String'));
h2 = findobj('Tag', 'editTolerance');
tolerance = str2double(get(h2,'String'));
h3 = findobj('Tag', 'checkboxShowplots');
showplots = get(h3,'Value');
h4 = findobj('Tag', 'chkCentre');
centred = get(h4,'Value');
types =[1];

hSX = findobj('Tag','editShiftx');
hSY = findobj('Tag','editShifty');
hFit = findobj('Tag','editFit');

scale = str2double(get(hScale,'String'));
shiftx = str2double(get(hSX,'String'));
shifty = str2double(get(hSY,'String'));
fit = str2double(get(hFit,'String'));
hC1 = findobj('Tag','editCell1');
cell1label = get(hC1, 'String');
hC2 = findobj('Tag','editCell2');
cell2label = get(hC2, 'String');
%Run analysis
hwb = waitbar(0,'Running analysis ...');
steps = 100;
step = 10;
waitbar(step / steps)
%show figs
axes(handles.axes2);
cla;
%Run analyser
N = NeuritesAnalyser(data.inputfile,roi,data.cell1,data.cell2);

title('Analysed ROI');
N = N.findSegments(minlength,tolerance);
step = step + 10;
waitbar(step / steps)

N = N.analyseSynapses(showplots,types);
step = step + 10;
waitbar(step / steps)

nocsvflag = 0;

if ( ~isempty(csvdata))
    csv1 = fullfile(csvdata.csvPath, csvdata.cell1file);
    csv2 = fullfile(csvdata.csvPath, csvdata.cell2file);
    N = N.measureSynapses(types,csv1,csv2,scale,shiftx,shifty,fit);
    pathname = csvdata.csvPath;
else
    nocsvflag = 1;
    pathname = data.imagePath;
end
step = step + 10;
waitbar(step / steps)

htable = findobj('Tag','uitableResults');
%Save data
[colnames,T] = N.generateTable(types,cell1label, cell2label);
set(htable,'data',T,'ColumnName',colnames);
%save to file
outputfile = fullfile(pathname, 'neurites_data.csv');
saveDataFile(outputfile, colnames, T);
%Save coords centred to DS centroid at [0,0]
if (centred)
    [colnames,T1] = N.generateCentredDataTable(types,cell1label, cell2label);
    %save to file
    outputfile = fullfile(pathname, 'neurites_data_centred.csv');
    saveDataFile(outputfile, colnames, T1);
end
%save analysis
dataN = struct('N', N, 'numsynapses',length(N.Synapses));
%hObject.UserData = dataN;
set(hObject,'UserData',dataN);
if (nocsvflag)
    csvmsg ='Load CSV files to get full statistics.';
else
    csvmsg = 'Statistics from CSV files.';
end
ctr = 0;
step = steps;
waitbar(step / steps)
close(hwb);

axes(handles.axes1);
hold on;
%Plot ROI for comparison
if (~isempty(roidata))
    plot(xi, yi, '--b','LineWidth', 1);
end
mycolors=['m' 'b' 'c'];
%mytypes =['En passant' 'End point' 'Intersection'];

for i=1:length(N.Synapses)
    syn = N.Synapses{i};
    if (ismember(syn.SynapseType,types))
        ctr = ctr + 1;
        %plot onto original
        plot(syn.MedianC1(:,1), syn.MedianC1(:,2),'color',...
            mycolors(syn.SynapseType),'marker','X','linestyle','none','LineWidth', 2);
        plot(syn.MedianC2(:,1), syn.MedianC2(:,2),'color',...
            mycolors(syn.SynapseType),'marker','X','linestyle','none','LineWidth', 2);
    end
end
hold off;

status = sprintf('Found %d synaptic regions. Saved to %s . %s', ctr,outputfile,csvmsg);
updateStatus(handles,status);
msgbox('Processing Complete!','Info');


function saveDataFile(outputfile,colnames,tabledata)
    T = array2table(tabledata);
    T.Properties.VariableNames = colnames(1,:);
    T
    writetable(T,outputfile);
    




% --- Executes on button press in btnSplitcells.
function btnSplitcells_Callback(hObject, eventdata, handles)
% hObject    handle to btnSplitcells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Use RGB image
h = findobj('Tag','btnBrowser');
%data = h.UserData;
data = get(h,'UserData');
I = data.img;
hBg = findobj('Tag','radioWhite');
bg = 256;
if (get(hBg,'Value') == 0)
    bg = 0;
end
[cell1,cell2] = SplitCells(I,bg);


% --- Executes on button press in btnReview.
function btnReview_Callback(hObject, eventdata, handles)
% hObject    handle to btnReview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = findobj('Tag','btnBrowser');
hData = get(h,'UserData');
I = hData.img;
f = figure('WindowStyle','normal');
im = imshow(I);
 
hold on;   
%Get Synapse data
hId = findobj('Tag','btnIdentify');
hNData = get(hId,'UserData');
types = [1];
hC1 = findobj('Tag','editCell1');
cell1label = strjoin(get(hC1, 'String'));
hC2 = findobj('Tag','editCell2');
cell2label = strjoin(get(hC2, 'String'));
c=2; %default cell2
if (~isempty(hNData) && hNData.numsynapses > 0)
    N1 = hNData.N;
%     changes = 0;
%     deleted = 0;
%     total = length(N1.Synapses);
    i = 1;
    syn = N1.Synapses{i};
    plot(N1.soma1.centroid(:,1), N1.soma1.centroid(:,2),'color', 'b',...
            'marker','o','linestyle','none','LineWidth', 2);
    plot(N1.soma2.centroid(:,1), N1.soma2.centroid(:,2),'color','b',...
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
       'Callback',{@setreviewtype,syn,N1});

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
       'Visible','off','Callback',{@acceptChanges,handles,cell1label, cell2label} );
   
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
   
  else
      msgbox('Cannot find synapse data! Run Analysis first.','Error');
end
function setcolour(source,callbackdata)
    c = source.Value;
    colours = ['b','r','y','c','m','k'];
    hId = findobj('Tag','btnReview');
    hData = get(hId,'UserData');
    hData.colour =colours(c);
    set(hId,'UserData',hData);
    
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
            p = plot(syn.SomapointsC1(:,1),syn.SomapointsC1(:,2),...
                'color',colour,'linestyle',linestyle,'LineWidth', 2);
        case 2
            p = plot(syn.EndpointsC1(:,1),syn.EndpointsC1(:,2),...
                'color',colour,'linestyle',linestyle,'LineWidth', 2);
            s1 = plot(syn.MedianC1(:,1), syn.MedianC1(:,2),'color',...
                'm','marker','X','linestyle','none','LineWidth', 2);
        case 3
            s1 = plot(syn.MedianC2(:,1), syn.MedianC2(:,2),'color',...
                'm','marker','X','linestyle','none','LineWidth', 2);
            p = plot(syn.SomapointsC2(:,1),syn.SomapointsC2(:,2), ...
                'color',colour,'linestyle',linestyle,'LineWidth', 2);
        case 4
            s1 = plot(syn.MedianC2(:,1), syn.MedianC2(:,2),'color',...
                'm','marker','X','linestyle','none','LineWidth', 2);
            p = plot(syn.EndpointsC2(:,1),syn.EndpointsC2(:,2), ...
                'color',colour,'linestyle',linestyle,'LineWidth', 2);
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
    
function changeSoma(source,callbackdata)
    nH = findobj('Tag', 'menu_File_loadimage');
    N1 = nH.UserData.analyser;
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
    
function clearplots(source,callbackdata)
    cla;
    h = findobj('Tag','btnBrowser');
    hData = get(h,'UserData');
    I = hData.img;
    %f = figure('WindowStyle','normal');
    im = imshow(I);
 
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

function acceptChanges(source,callbackdata,handles,cell1label, cell2label)
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
             updateStatus(handles,msg);
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
    
function editAnalysisfiles1_Callback(hObject, eventdata, handles)
% hObject    handle to editAnalysisfiles1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editAnalysisfiles1 as text
%        str2double(get(hObject,'String')) returns contents of editAnalysisfiles1 as a double


% --- Executes during object creation, after setting all properties.
function editAnalysisfiles1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editAnalysisfiles1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editTolerance_Callback(hObject, eventdata, handles)
% hObject    handle to editTolerance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTolerance as text
%        str2double(get(hObject,'String')) returns contents of editTolerance as a double


% --- Executes during object creation, after setting all properties.
function editTolerance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTolerance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMinlength_Callback(hObject, eventdata, handles)
% hObject    handle to editMinlength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMinlength as text
%        str2double(get(hObject,'String')) returns contents of editMinlength as a double


% --- Executes during object creation, after setting all properties.
function editMinlength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMinlength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxShowplots.
function checkboxShowplots_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxShowplots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxShowplots


% --- Executes on button press in checkboxEnpassant.
function checkboxEnpassant_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxEnpassant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxEnpassant


% --- Executes on button press in checkboxEndpoint.
function checkboxEndpoint_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxEndpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxEndpoint


% --- Executes on button press in checkboxIntersection.
function checkboxIntersection_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxIntersection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxIntersection




function editCell1_Callback(hObject, eventdata, handles)
% hObject    handle to editCell1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCell1 as text
%        str2double(get(hObject,'String')) returns contents of editCell1 as a double


% --- Executes during object creation, after setting all properties.
function editCell1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCell1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editCell2_Callback(hObject, eventdata, handles)
% hObject    handle to editCell2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCell2 as text
%        str2double(get(hObject,'String')) returns contents of editCell2 as a double


% --- Executes during object creation, after setting all properties.
function editCell2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCell2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
            




function editScale_Callback(hObject, eventdata, handles)
% hObject    handle to editScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editScale as text
%        str2double(get(hObject,'String')) returns contents of editScale as a double


% --- Executes during object creation, after setting all properties.
function editScale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to editScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editScale as text
%        str2double(get(hObject,'String')) returns contents of editScale as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to editScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editScale as text
%        str2double(get(hObject,'String')) returns contents of editScale as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editFit_Callback(hObject, eventdata, handles)
% hObject    handle to editScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editScale as text
%        str2double(get(hObject,'String')) returns contents of editScale as a double


% --- Executes during object creation, after setting all properties.
function editFit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editShiftx_Callback(hObject, eventdata, handles)
% hObject    handle to editScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editScale as text
%        str2double(get(hObject,'String')) returns contents of editScale as a double


% --- Executes during object creation, after setting all properties.
function editShiftx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editShifty_Callback(hObject, eventdata, handles)
% hObject    handle to editScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editScale as text
%        str2double(get(hObject,'String')) returns contents of editScale as a double


% --- Executes during object creation, after setting all properties.
function editShifty_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Menu_About_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_About (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Message = {'QBI Williams Neurites Annulus Analyser' ...
    'Description: Developed for the calculation of neurites area'...
    'under a stimulus region.' ...
    'Developed by Liz Cooper-Williams (e.cooperwilliams@uq.edu.au)'...
    '(c) Copyright QBI 2016' ...
    'Version: 1.x (Jun 2016)' ...
    'Source: https://github.com/QBI-Software/Neurites'};
msgbox(Message,'About')

% --------------------------------------------------------------------
function Menu_Help_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Message = {'Instructions for use' ...
    '1. Load image(tif or jpg)' ...
    '2. Load corresponding CSV data file' ...
    '3. Configure CSV overlay with Register CSV - adjust configuration values then Save config' ...
    '4. Add Annulus region (manual or auto)' ...
    '     a. Draw ROI with ROI tool and double click in region to save' ...
    '     b. Load from mask with File menu' ...
    '     c. Automatically apply with OD and ID from centroid'...
    '5. Run analysis' ...
    '6. Review results graphically - can adjust paths or delete synapses' ...
    '7. View output data in table which is also output to a csv file'};
h = msgbox(Message,'Instructions','Help');


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% --------------------------------------------------------------------
function Menu_SaveImage_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_SaveImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hB = findobj('Tag', 'menu_File_loadimage');
pathname = hB.UserData.csvPath;
if (isempty(pathname))
    pathname = '.';
end
saveas(handles.figure1, fullfile(pathname, 'neurites_img.png')); 




% --------------------------------------------------------------------
function Menu_LoadMask_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_LoadMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fileTypes = {  '*.tif;*.tiff;*.jpg;*.jpeg', 'Tiff/JPG images' };
[imageFile, imagePath, filterIndex] = ...
  uigetfile(fileTypes, ...
    'Select mask image', ...
    'MultiSelect', 'off');
inputPath = fullfile(imagePath,imageFile);
I = imread(inputPath);
[x,y] = size(I);
[r,c] = find(I > 0);
p = [r(1),c(1)]; %starting point
contour = bwtraceboundary(I,p,'E',8,Inf,'counterclockwise');
xi = contour(:,2);
yi = contour(:,1);
%save points
dataroi = struct('roi', I, 'xi',xi,'yi',yi,'x',x,'y',y);
hB = findobj('Tag', 'radioROI');
set(hB,'UserData',dataroi);
%display roi
axes(handles.axes2);
cla;
imshow(I)
hold on;
plot(xi,yi,'b','LineWidth',1);
hold off;
title('ROI Mask','FontSize', 10);
updateStatus(handles,'Mask loaded from file');
axes(handles.axes1);
hold on;
plot(xi,yi,'b--','LineWidth',1);
hold off;
%save ROI
hfiles = findobj('Tag','btnBrowser');
hfdata = get(hfiles,'UserData');
roifile = fullfile(hfdata.imagePath, 'neurites_roi.tif');
imwrite(I, roifile);


% --------------------------------------------------------------------
function Menu_LoadConfig_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_LoadConfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fileTypes = {  '*.csv', 'CSV files' };
[imageFile, imagePath, filterIndex] = ...
  uigetfile(fileTypes, ...
    'Select neurites config file', ...
    'MultiSelect', 'off');
configfile = fullfile(imagePath,imageFile);
if (exist(configfile, 'file') == 2)
    M = readtable(configfile);
    loadConfig(M,handles);
    msgbox('Config file loaded')
else
    msgbox('Unable to load config file')
end
    



function editOD_Callback(hObject, eventdata, handles)
% hObject    handle to editOD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editOD as text
%        str2double(get(hObject,'String')) returns contents of editOD as a double


% --- Executes during object creation, after setting all properties.
function editOD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editOD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editID_Callback(hObject, eventdata, handles)
% hObject    handle to editID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editID as text
%        str2double(get(hObject,'String')) returns contents of editID as a double


% --- Executes during object creation, after setting all properties.
function editID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function editCentroidX_Callback(hObject, eventdata, handles)
% hObject    handle to editCentroidX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCentroidX as text
%        str2double(get(hObject,'String')) returns contents of editCentroidX as a double


% --- Executes during object creation, after setting all properties.
function editCentroidX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCentroidX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editCentroidY_Callback(hObject, eventdata, handles)
% hObject    handle to editCentroidY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCentroidY as text
%        str2double(get(hObject,'String')) returns contents of editCentroidY as a double


% --- Executes during object creation, after setting all properties.
function editCentroidY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCentroidY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chkCentre.
function chkCentre_Callback(hObject, eventdata, handles)
% hObject    handle to chkCentre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    nH = findobj('Tag', 'menu_File_loadimage');
    N = nH.UserData.analyser;
    %Show coordinates
    hCX = findobj('Tag','editCentroidX');
    hCY = findobj('Tag','editCentroidY');
    if (isempty(get(hCX,'String')) || isempty(get(hCY,'String')))
        CentroidX = N.soma.centroid(1);
        CentroidY = N.soma.centroid(2);
    else
        CentroidX = str2double(get(hCX,'String'));
        CentroidY = str2double(get(hCY,'String'));
    end
    
    gca
    hold on
    p1 = plot(CentroidX, CentroidY,'color', 'b',...
                'marker','x','linestyle','none','LineWidth', 2);
    
    dcm_obj = datacursormode()
    set(dcm_obj,'DisplayStyle','datatip',...
    'SnapToDataVertex','off','Enable','on')
    status = sprintf('Click on figure at required centroid, then press Return.');
    updateStatus(handles,status);
    disp(status)
    
    % Wait while the user does this.
    pause 

    c_info = getCursorInfo(dcm_obj);
    pos = c_info.Position
    datacursormode off
    set(hCX, 'String',num2str(pos(1)));
	set(hCY, 'String',num2str(pos(2)));
    
    N.soma.centroid = [pos(1) pos(2)]
    
    delete(p1);
    plot(N.soma.centroid(:,1), N.soma.centroid(:,2),'color', 'r',...
            'marker','+','linestyle','none','LineWidth', 2);
    
    %save back
    nH.UserData.analyser = N;
    
    
