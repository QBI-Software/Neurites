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

% Last Modified by GUIDE v2.5 05-Jul-2017 11:13:48

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

function r = getProgramtype
    r = 'annulus';

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
%     hMagBox = immagbox(hFig,hIm);
%     pos = get(hMagBox,'Position');
%     set(hMagBox,'Position',[0 0 pos(3) pos(4)])
%     imoverview(hIm)
%     %Get the scroll panel API to programmatically control the view.
     api = iptgetapi(hSP);
%     %Get the current magnification and position.
%     %mag = api.getMagnification();
%     %r = api.getVisibleImageRect();
%     %View the top left corner of the image.
%     %api.setVisibleLocation(0.5,0.5)
%     %Change the magnification to the value that just fits.
     api.setMagnification(api.findFitMag())
%     %Zoom in to 1600% on the dark spot.
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
    updateStatus(handles,'Loading image, please wait ...');
    [imageFile, imagePath, ~] = ...
      uigetfile(fileTypes, ...
		'Select neuron image', ...
		'MultiSelect', 'off');
    inputPath = fullfile(imagePath,imageFile);
    
    I = imread(inputPath);
    %load image in panel
    axes(handles.axes1);
    hSP = get(handles.axes1,'parent');
    api = iptgetapi(hSP);
    api.replaceImage(I);
    % Reset magnification
    updateStatus(handles,'Loading image ....');
    mag = api.findFitMag()
    api.setMagnification(mag);
    htable = findobj('Tag','uitableResults');
    set(htable,'Data',[]);
    %Clear configuration
    clearConfig(getProgramtype, 0);
    clearCSVdata(getProgramtype);
    %clearROI( handles );
    %Save image analyser
    N = NeuritesAnnulusAnalyser(I); 
    %save image data
    data = struct('img',I, 'imagePath',imagePath,'inputfile', inputPath, 'analyser', N);
    set(hObject,'UserData',data);  
    %Notification
    updateStatus(handles,'Image loaded.');
    msgbox('Image loaded. Now load CSV file.')
    

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



% --------------------------------------------------------------------
function Menu_SaveImage_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_SaveImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hB = findobj('Tag', 'menu_File_loadimage');
pathname = hB.UserData.imagePath;
if (isempty(pathname))
    pathname = '.';
end
saveas(handles.figure1, fullfile(pathname, 'neurites_img.png')); 




% --------------------------------------------------------------------
function Menu_LoadAnnulusMask_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_LoadAnnulusMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fileTypes = {  '*.tif;*.tiff;*.jpg;*.jpeg', 'Tiff/JPG images' };
hB = findobj('Tag', 'menu_File_loadimage');
hfdata = get(hB,'UserData');
imagePath = hfdata.imagePath;
[imageFile, imagePath, ~] = ...
  uigetfile(fileTypes, ...
    'Select mask image', ...
    'MultiSelect', 'off',imagePath);
inputPath = fullfile(imagePath,imageFile);
I = imread(inputPath);
updateStatus(handles,'Annulus Mask loaded from file');
hfdata.roi_I=I;

%colors=['b' 'g' 'r' 'c' 'm' 'y'];
[B,~,~,~] = bwboundaries(I);
%display roi
axes(handles.axes1);
hold on;
for k=1:length(B),
  boundary = B{k};
  xi = boundary(:,2);
  yi = boundary(:,1);
  plot(xi,yi,'b--','LineWidth',1);
end

hold off;
hfdata.roi_xy = B
%save ROI as standard file (overwrites)
roifile = fullfile(hfdata.imagePath, 'neurites_annulus.tif');
if ~exist(roifile, 'file') == 2
    imwrite(I, roifile);
    hfdata.roifile = roifile;
end
%generate masked image
N = hfdata.analyser;
N = N.loadMask(I);
%imshow(N.maskedI);
hfdata.analyser = N;
hfdata.roi = 'annulus';
set(hB,'UserData',hfdata);

function Menu_LoadRectangleMask_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_LoadAnnulusMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fileTypes = {  '*.tif;*.tiff;*.jpg;*.jpeg', 'Tiff/JPG images' };
hB = findobj('Tag', 'menu_File_loadimage');
hfdata = get(hB,'UserData');
imagePath = hfdata.imagePath;
[imageFile, imagePath, ~] = ...
  uigetfile(fileTypes, ...
    'Select mask image', ...
    'MultiSelect', 'off',imagePath);
inputPath = fullfile(imagePath,imageFile);
I = imread(inputPath);
updateStatus(handles,'Rectangle Mask loaded from file');
hfdata.roi_I=I;

%colors=['b' 'g' 'r' 'c' 'm' 'y'];
[B,L,N,A] = bwboundaries(I);
%display roi
axes(handles.axes1);
hold on;
for k=1:length(B),
  boundary = B{k};
  xi = boundary(:,2);
  yi = boundary(:,1);
  plot(xi,yi,'b--','LineWidth',1);
end

hold off;
hfdata.roi_xy = B
%save ROI as standard file (overwrites)
roifile = fullfile(hfdata.imagePath, 'neurites_rectangle.tif');
if ~exist(roifile, 'file') == 2
    imwrite(I, roifile);
    hfdata.roifile = roifile;
end
%generate masked image
N = hfdata.analyser;
N = N.loadMask(I);
%imshow(N.maskedI);
hfdata.analyser = N;
hfdata.roi = 'rect';
set(hB,'UserData',hfdata);

% --------------------------------------------------------------------
function Menu_LoadConfig_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_LoadConfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fileTypes = {  '*.csv', 'CSV files' };
hB = findobj('Tag', 'menu_File_loadimage');
hfdata = get(hB,'UserData');
imagePath = hfdata.imagePath;
[imageFile, imagePath, ~] = ...
  uigetfile(fileTypes, ...
    'Select neurites config file', ...
    'MultiSelect', 'off',imagePath);
configfile = fullfile(imagePath,imageFile);
if (exist(configfile, 'file') == 2)
    M = readtable(configfile);
    loadConfig(M,getProgramtype,0);
    %set centroid
    N = hfdata.analyser;
    N.soma.centroid = [M.CentroidX M.CentroidY];
    hfdata.analyser = N;
    set(hB,'UserData',hfdata);
    msgbox('Config file loaded')
else
    msgbox('Unable to load config file')
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


       

% ---------------------RADIO BUTTONS--------------------------

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

% --- Executes on button press in radioAnnNone.
function radioAnnNone_Callback(hObject, eventdata, handles)
% hObject    handle to radioAnnNone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioAnnNone

% --- Executes on button press in radioDraw.
function radioDraw_Callback(hObject, eventdata, handles)
% hObject    handle to radioDraw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioDraw

% ---------------------PUSH BUTTONS--------------------------

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
    hW = findobj('Tag','editWidth');
    hH = findobj('Tag','editHeight');
    Scale = str2double(get(hScale,'String'));
    Fit = str2double(get(hFit,'String'));
    Shiftx = str2double(get(hSX,'String'));
    Shifty = str2double(get(hSY,'String'));
    CentroidX = str2double(get(hCx,'String'));
    CentroidY = str2double(get(hCy,'String'));
    AnnulusOD = str2double(get(hOD,'String'));
    AnnulusID = str2double(get(hID,'String'));
    Width= str2double(get(hW,'String'));
    Height = str2double(get(hH,'String'));

    T = table(Scale, Fit, Shiftx, Shifty, CentroidX, CentroidY, AnnulusOD, AnnulusID, Width, Height);
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
    
    
% --- Executes on button press in btnImtool.
function btnImtool_Callback(hObject, eventdata, handles)
% hObject    handle to btnImtool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = findobj('Tag','menu_File_loadimage');
%data = h.UserData;
data = get(h,'UserData');
imtool(data.img)
 

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
            [Scale, Shiftx, Shifty] = estimateOverlay(handles,I,csv1);
            CentroidX = N.soma.centroid(1);
            CentroidY = N.soma.centroid(2);
            AnnulusOD = 1;
            AnnulusID = 1;
            M = table(Scale, Shiftx, Shifty, CentroidX, CentroidY, AnnulusOD, AnnulusID);
        end
        loadConfig(M,getProgramtype,0);
        
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
          XC1(end+1) = (T1.EndX(i) * hscale) + shiftx;
          YC1(end+1) = (T1.EndY(i) * hscale) + shifty;  
          plot(XC1,-YC1,'color','c','LineStyle','-','LineWidth', 2);  
          XC1 = [];
          YC1 = [];  
        end
        
    end
    hold off;
end
             

% --- Executes on button press in btnAnnulus.
function btnAnnulus_Callback(hObject, eventdata, handles)
% hObject    handle to btnAnnulus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of btnAnnulus
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
        mask = applyAnnulus(id, od, scale, centroid);
        N = N.loadMask(mask);
        imshow(N.maskedI);
        roifile = fullfile(data.imagePath, 'neurites_annulus.tif');
        imwrite(mask, roifile);
        status = sprintf('Annulus mask created: %s.', roifile);
        %Update analyser data
        data.analyser = N;
        data.roi = 'annulus';
        set(h,'UserData',data);
        updateStatus(handles,status);
        disp(status)
        
    end
end


% --- Executes on button press in btnAnalysis.
function btnAnalysis_Callback(hObject, eventdata, handles)
    h = findobj('Tag','menu_File_loadimage');
    %data = h.UserData;
    data = get(h,'UserData');
    N = data.analyser;
    hCSV = findobj('Tag','Menu_File_loadcsv');
    csvdata = get(hCSV,'UserData');
    updateStatus(handles,'Running Analysis ...');
    if ( ~isempty(csvdata) && ~isempty(N.maskedI))
        csv = fullfile(csvdata.csvPath, csvdata.csvFile);
        CSVFile = readtable(csv);
        %Get scale factor
        hScale = findobj('Tag','editScale');
        Scale = str2num(get(hScale,'String'));
        hX = findobj('Tag','editShiftx');
        Shiftx = str2num(get(hX,'String'));
        hY = findobj('Tag','editShifty');
        Shifty = str2num(get(hY,'String'));
        hW = findobj('Tag','choiceDirection');
        d1 = get(hW,'UserData');
        direction = lower(d1(get(hW,'Value')));
        direction = direction{1}; %string from cellstr
        hID = findobj('Tag','editID');
        id = str2double(get(hID,'String'));
        hOD = findobj('Tag','editOD');
        od = str2double(get(hOD,'String'));
        hAL = findobj('Tag','editLength');
        arclength = str2double(get(hAL,'String'));
        if (isempty(arclength)||isnan(arclength))
            arclength=45;
        end
        hML = findobj('Tag','editMidline');
        midline = str2double(get(hML,'String'));
        if (isempty(midline)||isnan(midline))
            midline=45;
        end
        hSr = findobj('Tag','chkSingleRun');
        if (get(hSr,'Value') > 0)
            single = 1;      
        else
            single = 0;
        end
        if (~isempty(data.roi) && strcmp(data.roi,'annulus')>0)
            [status,regionMap] = analyseAnnulus(Scale,Shiftx,Shifty,N, (od * Scale)/2,(id * Scale)/2, midline, arclength, single, CSVFile);
        elseif (~isempty(data.roi) && strcmp(data.roi,'rect')>0)
            [status,regionMap] = analyseRectangle(Scale,Shiftx,Shifty,N, direction, single, CSVFile);
        else
            msgbox('Please create Annulus or Rectangle mask first')
            return
        end
        N.Regions = regionMap;
        %Generate table
        arcs = keys(regionMap)
        k=0;
        for i=1:numel(arcs)
            n = regionMap(arcs{i});
            if (length(arcs)==1)
                n = n{1,1};
            end
            for j = 1:numel(n.neurites)
                k = k+1
                neurite = n.neurites(j)
                lineStruct(k,1).ROI = n.id;                         %ID of ROI
                lineStruct(k,1).ROIPosition = n.midline;            %top left x position OR midline of arc in degrees eg 45o
                lineStruct(k,1).ROILength = n.slength;              %length of arc in degrees eg 45o
                lineStruct(k,1).ROIArea = n.sarea;                  %area of arc (um2) 
                lineStruct(k,1).CSVRow = n.neurites(j).crow;        %corresponding row number in CSV
                
                lineStruct(k,1).Tree = n.neurites(j).tree;          %tree number of this neurite
                if (length(n.neurites(j).branch) == 1)
                    lineStruct(k,1).Branch = n.neurites(j).branch;  %first branch number/s of this neurite
                else
                    lineStruct(k,1).Branch = n.neurites(j).branch(1);
                end
                lineStruct(k,1).BranchLength = n.neurites(j).blength;
                lineStruct(k,1).Length = n.neurites(j).nlength;     %length from min to max points (should be total neurite segment length) (um)                
                lineStruct(k,1).X = n.neurites(j).xy(1);
                lineStruct(k,1).Y = n.neurites(j).xy(2);
                lineStruct(k,1).SomaDist = n.neurites(j).somad;     %distance back to soma (um) from x,y
                %lineStruct(k,1).SomaVol = n.neurites(j).somav;      %Volume of dendrite back to soma (um3) ie maxpoint
                %lineStruct(k,1).SomaSA = n.neurites(j).somas;       %SA of dendrite back to soma (um2) ie maxpoint
                %lineStruct(k,1).Color = n.color; 
            end
        end
    
        T = struct2table(lineStruct(1:end,:));
        colnames = T.Properties.VariableNames
        htable = findobj('Tag','uitableResults');
        % Data must be a numeric, logical, or cell array
        %colnames = {'ArcMidline' 'Area' 'Tree' 'Branch' 'Length' 'NArea' 'Somax'};
        set(htable,'data',[T.ROI,T.ROIPosition,T.ROILength,T.ROIArea,T.CSVRow,T.Tree,T.Branch,T.BranchLength,T.Length,T.X,T.Y,T.SomaDist],'ColumnName',colnames);
        
        pathname = data.imagePath;
        %save to file
        outputfile = fullfile(pathname, 'neurites_annulus_data.csv');
        if single
            %append to table data
            
            if exist(outputfile, 'file') == 2
                T1 = readtable(outputfile);
                Tnew = [T1;T];
            else
                Tnew = T;
            end
            writetable(Tnew,outputfile);
        else
            writetable(T,outputfile);
        end
        status = sprintf('Analysis complete. %s. Table written to %s',status, outputfile);
        updateStatus(handles,status);
    else
        msgbox('Create Annulus mask first');
    end
    
% --- Executes on button press in chkSingleRun.
function chkSingleRun_Callback(hObject, eventdata, handles)
% hObject    handle to chkSingleRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnChangeCentroid.
function btnChangeCentroid_Callback(hObject, eventdata, handles)
% hObject    handle to btnChangeCentroid (see GCBO)
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
    
    plotCentroid(CentroidX, CentroidY);
    
    dcm_obj = datacursormode()
    set(dcm_obj,'DisplayStyle','datatip',...
    'SnapToDataVertex','off','Enable','on')
    status = sprintf('Click on figure at required centroid, then press Return.');
    msgbox(status)
    
    % Wait while the user does this.
    pause 

    c_info = getCursorInfo(dcm_obj);
    pos = c_info.Position
    datacursormode off
    set(hCX, 'String',num2str(pos(1)));
	set(hCY, 'String',num2str(pos(2)));
    
    N.soma.centroid = [pos(1) pos(2)]
    p1 = plotCentroid(N.soma.centroid(:,1), N.soma.centroid(:,2));
    
    %delete(p1);
    nC = findobj('Tag', 'btnCentroid_Callback');
    
    if ~isempty(get(nC,'UserData'))
        p0 = nC.UserData.centroid;
        delete(p0);
        nC.UserData.centroid = p1;
    else
        data = struct('centroid',p1);
        set(nC,'UserData',data);
    end
       
    %save back
    nH.UserData.analyser = N;
    set(nH, 'UserData', nH.UserData);

    
% --- Executes on button press in btnCentroid.
function btnCentroid_Callback(hObject, eventdata, handles)
% hObject    handle to btnCentroid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of btnCentroid
nH = findobj('Tag', 'menu_File_loadimage');
N = nH.UserData.analyser;
CentroidX = N.soma.centroid(1);
CentroidY = N.soma.centroid(2);
pc = plotCentroid(CentroidX, CentroidY); 
data = struct('centroid',pc);
set(hObject,'UserData',data); 


function saveDataFile(outputfile,colnames,tabledata)
    T = array2table(tabledata);
    T.Properties.VariableNames = colnames(1,:);
    T
    writetable(T,outputfile);
    
% ---------------------TEXT FIELDS--------------------------

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



% --- Executes on slider movement.
function sliderMidline_Callback(hObject, eventdata, handles)
% hObject    handle to sliderMidline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

nH = findobj('Tag', 'editMidline');
set(nH, 'String', get(hObject,'Value'));

% --- Executes during object creation, after setting all properties.
function sliderMidline_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderMidline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderLength_Callback(hObject, eventdata, handles)
% hObject    handle to sliderLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
nH = findobj('Tag', 'editLength');
set(nH, 'String', get(hObject,'Value'));

% --- Executes during object creation, after setting all properties.
function sliderLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function editMidline_Callback(hObject, eventdata, handles)
% hObject    handle to editMidline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMidline as text
%        str2double(get(hObject,'String')) returns contents of editMidline as a double
nH = findobj('Tag', 'sliderMidline');
set(nH, 'Value', str2double(get(hObject,'String')));

% --- Executes during object creation, after setting all properties.
function editMidline_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMidline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editLength_Callback(hObject, eventdata, handles)
% hObject    handle to editLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editLength as text
%        str2double(get(hObject,'String')) returns contents of editLength as a double
nH = findobj('Tag', 'sliderLength');
set(nH, 'Value', str2double(get(hObject,'String')));

% --- Executes during object creation, after setting all properties.
function editLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in btnReloadImage.
function btnReloadImage_Callback(hObject, eventdata, handles)
% hObject    handle to btnReloadImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    axes(handles.axes1);
    hfiles = findobj('Tag','menu_File_loadimage');
    hfdata = get(hfiles,'UserData');
    I = hfdata.img;
    cla;
    imshow(I);
    updateStatus(handles,'Reloaded image');


% --- Executes on button press in pbRectangle.
function pbRectangle_Callback(hObject, eventdata, handles)
% hObject    handle to pbRectangle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    axes(handles.axes1);
    hfiles = findobj('Tag','menu_File_loadimage');
    hfdata = get(hfiles,'UserData');
    N = hfdata.analyser;
    hScale = findobj('Tag','editScale');
    scale = str2double(get(hScale,'String'));
    hfRectWidth = findobj('Tag','editHeight');
    hfRectHeight = findobj('Tag','editWidth');
    height = str2double(get(hfRectHeight,'String'));
    width = str2double(get(hfRectWidth,'String'));
    centroid =N.soma.centroid;
    mask = applyRectangle(width,height, scale, centroid);
    N = N.loadMask(mask);
    imshow(N.maskedI);
    roifile = fullfile(hfdata.imagePath, 'neurites_rectangle.tif');
    imwrite(mask, roifile);
    status = sprintf('Rectangle mask created: %s.', roifile);
    %Update analyser data
    hfdata.analyser = N;
    set(hfiles,'UserData',hfdata);
    updateStatus(handles,status);
    




function editWidth_Callback(hObject, eventdata, handles)
% hObject    handle to editWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editWidth as text
%        str2double(get(hObject,'String')) returns contents of editWidth as a double


% --- Executes during object creation, after setting all properties.
function editWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editHeight_Callback(hObject, eventdata, handles)
% hObject    handle to editHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editHeight as text
%        str2double(get(hObject,'String')) returns contents of editHeight as a double


% --- Executes during object creation, after setting all properties.
function editHeight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editFit_Callback(hObject, eventdata, handles)
% hObject    handle to editFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFit as text
%        str2double(get(hObject,'String')) returns contents of editFit as a double


% --- Executes during object creation, after setting all properties.
function editFit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rbAnn.
function rbAnn_Callback(hObject, eventdata, handles)
% hObject    handle to rbAnn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbAnn
hfiles = findobj('Tag','menu_File_loadimage');
hfdata = get(hfiles,'UserData');
if (get(hObject,'Value') > 0)
    hfdata.roi = 'annulus';
else
    hfdata.roi = 'rect';
end
set(hfiles,'UserData', hfdata);



% --- Executes on button press in rbRect.
function rbRect_Callback(hObject, eventdata, handles)
% hObject    handle to rbRect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbRect
hfiles = findobj('Tag','menu_File_loadimage');
hfdata = get(hfiles,'UserData');
if (get(hObject,'Value') > 0)
    hfdata.roi = 'rect';
else
    hfdata.roi = 'annulus';
end
set(hfiles,'UserData', hfdata);



% --- Executes on selection change in choiceDirection.
function choiceDirection_Callback(hObject, eventdata, handles)
% hObject    handle to choiceDirection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns choiceDirection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from choiceDirection


% --- Executes during object creation, after setting all properties.
function choiceDirection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to choiceDirection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
   
    
end



% --- Executes on button press in btnCSVAnalysis.
function btnCSVAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to btnCSVAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    h = findobj('Tag','menu_File_loadimage');
    %data = h.UserData;
    data = get(h,'UserData');
     
    hCSV = findobj('Tag','Menu_File_loadcsv');
    csvdata = get(hCSV,'UserData');
    updateStatus(handles,'Running Analysis ...');
    if ( ~isempty(csvdata))
        csv = fullfile(csvdata.csvPath, csvdata.csvFile);
        
        %Get scale factor
        hScale = findobj('Tag','editScale');
        scale = 1;%str2num(get(hScale,'String'));
        hA = findobj('Tag','rbAnn');
        annulus = get(hA,'Value');
       
        if (annulus)
            %annulus
            hID = findobj('Tag','editID');
            id = str2double(get(hID,'String'));
            hOD = findobj('Tag','editOD');
            od = str2double(get(hOD,'String'));
            shape = [id, od];
            annulus = true;
        else
            %rectangle
            annulus = false;
            hID = findobj('Tag','editWidth');
            width = str2double(get(hID,'String'));
            hID = findobj('Tag','editHeight');
            height = str2double(get(hID,'String'));
            shape = [width,height];
            hW = findobj('Tag','choiceDirection');
            d1 = get(hW,'UserData');
            direction = lower(d1(get(hW,'Value')));
            direction = direction{1}; %string from cellstr
            hSr = findobj('Tag','chkSingleRun');
            if (get(hSr,'Value') > 0)
                single = 1;      
            else
                single = 0;
            end
        end
        
        N = NeuritesCSVAnnulusAnalyser(csv, annulus, shape, scale);
        updateStatus(handles,'Total regions found=%d', length(N.regions));
        [colnames,tabledata] = N.generateTable();
        htable = findobj('Tag','uitableResults');
        set(htable,'data', tabledata, 'ColumnName', colnames);
        outputfile = fullfile(csvdata.csvPath, 'neurites_csv_annulus.csv');
        saveDataFile(outputfile,colnames,tabledata);
    else
        updateStatus(handles,'ERROR: CSV file not loaded');
    end

        
        

