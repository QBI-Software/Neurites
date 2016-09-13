function varargout = NeuritesAppUI(varargin)
% NEURITESAPPUI MATLAB code for NeuritesAppUI.fig
%      NEURITESAPPUI, by itself, creates a new NEURITESAPPUI or raises the existing
%      singleton*.
%
%      H = NEURITESAPPUI returns the handle to a new NEURITESAPPUI or the handle to
%      the existing singleton*.
%
%      NEURITESAPPUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NEURITESAPPUI.M with the given input arguments.
%
%      NEURITESAPPUI('Property','Value',...) creates a new NEURITESAPPUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NeuritesAppUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NeuritesAppUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% ***NB***% For R2014a and earlier: 
%       data = get(h,'UserData');
%       set(hObject,'UserData',data); 
% Edit the above text to modify the response to help NeuritesAppUI

% Last Modified by GUIDE v2.5 15-Jan-2016 11:28:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NeuritesAppUI_OpeningFcn, ...
                   'gui_OutputFcn',  @NeuritesAppUI_OutputFcn, ...
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
    r = 'synapse';

    
% --- Executes just before NeuritesAppUI is made visible.
function NeuritesAppUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NeuritesAppUI (see VARARGIN)

% Choose default command line output for NeuritesAppUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using NeuritesAppUI.
if strcmp(get(hObject,'Visible'),'off')
    %(rand(5));
    %load image in panel
    axes(handles.axes1);
    hFig = figure(get(handles.axes1,'parent'));%handles.axes1.Parent;
    hIm = imshow('synapse_EM.jpg');
    hSP = imscrollpanel(hFig,hIm); % Handle to scroll panel.
    set(hSP,'Units','pixels',...
        'Position',[50 289 551 351]) %match with GUI box
	    
    % Add a Magnification Box and an Overview tool.
    hMagBox = immagbox(hFig,hIm);
    pos = get(hMagBox,'Position');
    set(hMagBox,'Position',[0 0 pos(3) pos(4)])
    imoverview(hIm)
    %Get the scroll panel API to programmatically control the view.
    api = iptgetapi(hSP);
    %Get the current magnification and position.
    mag = api.getMagnification();
    r = api.getVisibleImageRect();
    %View the top left corner of the image.
    %api.setVisibleLocation(0.5,0.5)
    %Change the magnification to the value that just fits.
    api.setMagnification(api.findFitMag())
    %Zoom in to 1600% on the dark spot.
    api.setMagnificationAndCenter(16,306,800)
    %test
    I = imread('synapse_EM.jpg');
    api.replaceImage(I)
    
end

% UIWAIT makes NeuritesAppUI wait for user response (see UIRESUME)
 %uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = NeuritesAppUI_OutputFcn(hObject, eventdata, handles)
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

% --------------------------------------------------------------------
function Menu_About_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_About (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Message = {'QBI Williams Neurites Synapse Analyser' ...
    'Description: Developed for Prof Stephen Williams and Dr Simon de Croft (Williams lab) '...
    'for the detection and analysis of synaptic regions of two neurons.' ...
    'Developed by Liz Cooper-Williams (e.cooperwilliams@uq.edu.au)'...
    '(c) Copyright QBI 2015' ...
    'Version: 1.x (Sep 2015)' ...
    'Version: 2.x (Oct 2015)' ...
    'Source: https://github.com/QBI-Software/Neurites'};
msgbox(Message,'About')

% --------------------------------------------------------------------
function Menu_Help_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Message = {'Instructions for use' ...
    '1. Load image(tif or jpg), check histogram and split cells' ...
    '2. Load 2 CSV data files for each cell and check labels (eg DS.csv and SBAC.csv)' ...
    '3. Configure CSV overlay with Register CSV - adjust configuration values then Save config' ...
    '4. Draw ROI with ROI tool and double click in region to save (or load from mask with File menu)' ...
    '5. Run analysis' ...
    '6. Review results graphically - can adjust paths or delete synapses' ...
    '7. Compass - shows plot of angles between synapses and cell somas' ...
    '8. View output data in table and csv file'};
h = msgbox(Message,'Instructions','Help');


% --------------------------------------------------------------------
function menu_File_openImg_Callback(hObject, eventdata, handles)
% hObject    handle to menu_File_openImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hObject2 = findobj('Tag', 'btnBrowser');
NeuritesAppUI('btnBrowser_Callback',hObject2,eventdata,handles)


% --------------------------------------------------------------------
function Menu_SaveImage_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_SaveImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hB = findobj('Tag', 'btnBrowser');
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
    loadConfig(M,getProgramtype,0);
    msgbox('CSV files loaded - now run Register CSV and adjust config')
else
    msgbox('Unable to load config file')
end
    


% --------------------------------------------------------------------
function Menu_File_Opencsv_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_File_Opencsv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hObject2 = findobj('Tag', 'btnAnalysisFiles');
NeuritesAppUI('btnAnalysisFiles_Callback',hObject2,eventdata,handles)



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



function txtInputfilename_Callback(hObject, eventdata, handles)
% hObject    handle to txtInputfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtInputfilename as text
%        str2double(get(hObject,'String')) returns contents of txtInputfilename as a double


% --- Executes during object creation, after setting all properties.
function txtInputfilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtInputfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnBrowser.
function btnBrowser_Callback(hObject, eventdata, handles)
% hObject    handle to btnBrowser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    fileTypes = {  '*.tif;*.tiff;*.jpg;*.jpeg', 'Tiff/JPG images' };
    
    [imageFile, imagePath, ~] = ...
      uigetfile(fileTypes, ...
		'Select neuron image', ...
		'MultiSelect', 'off');
    inputPath = fullfile(imagePath,imageFile);
    set(handles.txtInputfilename,'string',inputPath)
    I = imread(inputPath);
    %ig = rgb2gray(I);
    %load image in panel
    axes(handles.axes1);
    %hFig = figure(get(handles.axes1,'parent'));%handles.axes1.Parent;
    hSP = get(handles.axes1,'parent');
    %hSP = handles.axes1.Parent; %scrollpanel
    %hFig = hSP.Parent; %figure
    api = iptgetapi(hSP);
    api.replaceImage(I)
    % Reset magnification
    mag = api.findFitMag()
    api.setMagnification(mag);
    %Panel 2 - show histogram  
    hBg = findobj('Tag','radioWhite');
    bg = 256;
    if (get(hBg,'Value') == 0)
        bg = 0;
    end
    [cell1,cell2] = SplitCells(I,bg);
    data = struct('img',I, 'imagePath',imagePath,'inputfile', inputPath, ...
        'cell1',cell1,'cell2',cell2);
    set(hObject,'UserData',data); 
    %hObject.UserData = data;
    %clear any previous data
    htable = findobj('Tag','uitableResults');
    set(htable,'Data',[]);
    %or set(htable,'Visible','off')
    %Clear configuration
    clearConfig(getProgramtype,0);
    clearCSVdata(getProgramtype);
    %Directions
    updateStatus(handles,'Image loaded. Check split cells correspond to Cell1 and Cell2 labels');


% --- Executes on button press in btnIdentify.
function btnIdentify_Callback(hObject, eventdata, handles)
% hObject    handle to btnIdentify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%clearplots();
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
%cla;
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
    %Get Synapse data
    hId = findobj('Tag','btnIdentify');
    hNData = get(hId,'UserData');
    N = hNData.N;
    types = [1];
    hC1 = findobj('Tag','editCell1');
    cell1label = strjoin(get(hC1, 'String'));
    hC2 = findobj('Tag','editCell2');
    cell2label = strjoin(get(hC2, 'String'));
    if (~isempty(hNData) && hNData.numsynapses > 0)
        NeuritesReview(I,N,cell1label,cell2label);
    else
        msgbox('Cannot find synapse data! Run Analysis first.','Error');
    end
    
    
    
% --- Executes on button press in btnAnalysisFiles.
function btnAnalysisFiles_Callback(hObject, eventdata, handles)
% hObject    handle to btnAnalysisFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    fileTypes = {  '*.csv', 'CSV' };
    hCf1 = findobj('Tag', 'editCell1');
    cell1ID=get(hCf1,'String');% 'DS';
    if (iscell(cell1ID))
        cell1ID = strjoin(cell1ID);
    end
    
    hfiles = findobj('Tag','btnBrowser');
    %hfdata = hfiles.UserData;
    hfdata = get(hfiles,'UserData');
    csvPath = hfdata.imagePath;
    [csvFile, csvPath, ~] = ...
      uigetfile(fileTypes, ...
		'Select analysis files', ...
		'MultiSelect', 'off',csvPath);
    cell1 = csvFile;
    set(handles.editAnalysisfiles1,'string',csvFile);
    csv1 = fullfile(csvPath, csvFile);
    status = sprintf('Cell 1: %s', csv1);
    updateStatus(handles,status);
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
        data = struct('csvPath', csvPath, 'cell1file', cell1, 'cell1ID',cell1ID);
    else
        data.csvPath = csvPath;
        data.cell1file = cell1;
        data.cell1ID = cell1ID;
    end
    %hObject.UserData = data;
    set(hObject,'UserData',data);
    if (isfield(data,'cell1file') && isfield(data,'cell2file'))
         msgbox('CSV files loaded. Run configuration with Register CSV.')
    end
    
    % --- Executes on button press in btnAnalysisFiles2.
function btnAnalysisFiles2_Callback(hObject, eventdata, handles)
% hObject    handle to btnAnalysisFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    fileTypes = {  '*.csv', 'CSV' };
    hCf2 = findobj('Tag', 'editCell2');
    cell2ID=get(hCf2,'String');% 'SBAC';
    if (iscell(cell2ID))
        cell2ID = strjoin(cell2ID);
    end
    hcsv1 = findobj('Tag','btnAnalysisFiles');
    hfiles = findobj('Tag','btnBrowser');
    %hfdata = hfiles.UserData;
    hfdata = get(hfiles,'UserData');
    csvPath = hfdata.imagePath;
    [csvFile, csvPath, ~] = ...
      uigetfile(fileTypes, ...
		'Select analysis files', ...
		'MultiSelect', 'off',csvPath);
    cell2 = csvFile;
    set(handles.editAnalysisfiles2,'string',csvFile);
    csv2 = fullfile(csvPath, csvFile);
    status = sprintf('Cell 2: %s', csv2);
    updateStatus(handles,status);
    %Check CSV headers
    hdrs = 'Tree,Order,StartX,StartY,StartZ,EndX,EndY,EndZ,PointType,Length(µm),LengthToBeginning(µm)';
    r1 = checkHeaders(csv2,hdrs);
    if (length(r1) > 1)
        errhdr = sprintf('ERROR: CSV file headers need to match: %s. \nFields incorrect:%s %s', hdrs,r1);
        errordlg(errhdr,'CSV Files error')
        return 
    end
    %save to userdata
    data = get(hObject,'UserData');
    if (isempty(data))
        data = struct('csvPath', csvPath, 'cell2file', cell2, 'cell2ID',cell2ID);
    else
        data.csvPath = csvPath;
        data.cell2file = cell2;
        data.cell2ID = cell2ID;
    end
    %hObject.UserData = data;
    set(hcsv1,'UserData',data);
    if (isfield(data,'cell1file') && isfield(data,'cell2file'))
         msgbox('CSV files loaded. Run configuration with Register CSV.')
    end
    
% --- Executes on button press in btnRegister.
function btnRegister_Callback(hObject, eventdata, handles)
% hObject    handle to btnRegister (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = findobj('Tag','btnBrowser');
%data = h.UserData;
data = get(h,'UserData');
I = data.cell1;
%csvdata = handles.btnAnalysisFiles.UserData;
hCSV = findobj('Tag','btnAnalysisFiles');
csvdata = get(hCSV,'UserData');

if ( ~isempty(csvdata))
    csv1 = fullfile(csvdata.csvPath, csvdata.cell1file);
    %Read values
    hScale = findobj('Tag','editScale');
    hSX = findobj('Tag','editShiftx');
    hSY = findobj('Tag','editShifty');
    Fit = 10;
    Minlength = 10;
    Tolerance = 1;
    if isempty(get(hScale,'String'))
        %load config data if empty
        configfile =fullfile(csvdata.csvPath, 'neurites_config.csv');
        if (exist(configfile, 'file') == 2)
            M = readtable(configfile);
        else
            [Scale, Shiftx, Shifty] = estimateOverlay(handles,I, csv1);
            Cell1 = csvdata.cell1ID;
            Cell2 = csvdata.cell2ID;
            M = table(Cell1,Cell2, Scale, Fit, Shiftx, Shifty, Minlength,Tolerance);
           
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
    axes(handles.axes2);
    imshow(I);
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


% --- Executes on button press in btnCompass. %TODO
function btnCompass_Callback(hObject, eventdata, handles)
% hObject    handle to btnCompass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 hId = findobj('Tag','btnIdentify');
 hNData = get(hId,'UserData');
 if (isempty(hNData))
     msgbox('No data found - please run Analysis first');
 else
     N = hNData.N;
     l = length(N.Synapses);
     theta1 = [l];
     theta2 = [l];
     rho1 = [l];
     rho2 = [l];
     %display roi
     figure
     %hold on;
     for i=1:l
         syn = N.Synapses{i};
         theta1(end+1) = syn.ThetaC1;
         theta2(end+1) = syn.ThetaC2;
         rho1(end+1) = syn.RhoC1;
         rho2(end+1) = syn.RhoC2;
         
     end
     
     [x1,y1] = pol2cart(2 * pi - theta1,rho1);
     [x2,y2] = pol2cart(2 * pi - theta2,rho2);
     
     hCell1 = findobj('Tag','editCell1');
     hCell2 = findobj('Tag','editCell2');
     %Add legend relative to cell label
     if (strcmp(get(hCell1, 'String'),'DS'))
        color1 = '-r';
        color2 = '-g';
     else
        color1 = '-g';
        color2 = '-r';
     end
     
     compass(x1,-y1,color1);
     hold on
     compass(x2,-y2,color2);
     set(gca,'View',[-90 90],'YDir','reverse');
     legend([get(hCell1, 'String'),get(hCell2, 'String')]);
 end


% --- Executes on button press in btnRose.
function btnRose_Callback(hObject, eventdata, handles)
% hObject    handle to btnRose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 hId = findobj('Tag','btnIdentify');
 hNData = get(hId,'UserData');
 if (isempty(hNData))
     msgbox('No data found - please run Analysis first');
 else
     N = hNData.N;
     l = length(N.Synapses);
     theta1 = [l];
     theta2 = [l];
     rho1 = [l];
     rho2 = [l];
     %display roi
     figure
     %hold on;
     for i=1:l
         syn = N.Synapses{i};
         theta1(end+1) = syn.ThetaC1;
         theta2(end+1) = syn.ThetaC2;
         rho1(end+1) = syn.RhoC1;
         rho2(end+1) = syn.RhoC2;
         
     end
     
     hCell1 = findobj('Tag','editCell1');
     hCell2 = findobj('Tag','editCell2');
          
     %show plots
     
     rose(theta1);
     hold on
     rose(theta2);
     set(gca,'View',[-90 90],'YDir','reverse');
     legend([get(hCell1, 'String'),get(hCell2, 'String')])
 end
    

% --- Executes on button press in saveConfig.
function saveConfig_Callback(hObject, eventdata, handles)
% hObject    handle to saveConfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    %Results table
    hCell1 = findobj('Tag','editCell1');
    hCell2 = findobj('Tag','editCell2');
    hScale = findobj('Tag','editScale');
    hSX = findobj('Tag','editShiftx');
    hSY = findobj('Tag','editShifty');
    hFit = findobj('Tag','editFit');
    hMin = findobj('Tag','editMinlength');
    hTol = findobj('Tag','editTolerance');
    Cell1 = get(hCell1,'String');
    Cell2 = get(hCell2,'String');
    Scale = str2double(get(hScale,'String'));
    Fit = str2double(get(hFit,'String'));
    Shiftx = str2double(get(hSX,'String'));
    Shifty = str2double(get(hSY,'String'));
    Minlength = str2double(get(hMin,'String'));
    Tolerance = str2double(get(hTol,'String'));

    T = table(Cell1, Cell2, Scale, Fit, Shiftx, Shifty, Minlength,Tolerance);
    %save to file
    %csvdata = handles.btnAnalysisFiles.UserData;
    hCSV = findobj('Tag','btnAnalysisFiles');
    csvdata = get(hCSV,'UserData');
    outputfile = 'neurites_config.csv';
    if ( ~isempty(csvdata))
        outputfile =fullfile(csvdata.csvPath, outputfile);
    end
    writetable(T,outputfile);
    status = sprintf('Config file saved:%s',outputfile);
    updateStatus(handles,status);
    
     

% --------------------------------------------------------------------
% --- Executes on button press in btnImtool.
function btnImtool_Callback(hObject, eventdata, handles)
% hObject    handle to btnImtool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = findobj('Tag','btnBrowser');
%data = h.UserData;
data = get(h,'UserData');
imtool(data.img)


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
    %save data
    dataroi = struct('roi', BW, 'xi',xi,'yi',yi,'x',x,'y',y);
    hfiles = findobj('Tag','btnBrowser');
    %hfdata = hfiles.UserData;
    hfdata = get(hfiles,'UserData');
    roifile = fullfile(hfdata.imagePath, 'neurites_roi.tif');
    imwrite(BW, roifile);
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
    [RCols, RData] = areaPolyROI(dataroi,scale);
    %Results table
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


function editAnalysisfiles2_Callback(hObject, eventdata, handles)
% hObject    handle to editAnalysisfiles2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editAnalysisfiles2 as text
%        str2double(get(hObject,'String')) returns contents of editAnalysisfiles2 as a double


% --- Executes during object creation, after setting all properties.
function editAnalysisfiles2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editAnalysisfiles2 (see GCBO)
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

function clearplots(source,callbackdata)
    cla;
    h = findobj('Tag','btnBrowser');
    hData = get(h,'UserData');
    I = hData.img;
    %f = figure('WindowStyle','normal');
    im = imshow(I);

 
