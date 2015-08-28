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

% Edit the above text to modify the response to help NeuritesAppUI

% Last Modified by GUIDE v2.5 20-Aug-2015 16:08:53

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
    hFig = handles.axes1.Parent;
    hIm = imshow('tmp_roi.tif');
    hSP = imscrollpanel(hFig,hIm); % Handle to scroll panel.
    set(hSP,'Units','normalized',...
        'Position',[0.035 .254 0.393 .597]) %match with GUI box
	%title('Default');
    %setTitlePlot1('Default');
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
    I = imread('example.tif');
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

% --- Executes on button press in btnUpdate.
function btnUpdate_Callback(hObject, eventdata, handles)
% hObject    handle to btnUpdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
cla;

popup_sel_index = get(handles.popupmenu1, 'Value');
switch popup_sel_index
    case 1
        plot(rand(5));
    case 2
        plot(sin(1:0.01:25.99));
    case 3
        bar(1:.5:10);
    case 4
        plot(membrane);
    case 5
        surf(peaks);
end


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

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

    fileTypes = {  '*.tif;*.tiff', 'Tiff-image' };
    
    [imageFile, imagePath, filterIndex] = ...
      uigetfile(fileTypes, ...
		'Select neuron image', ...
		'MultiSelect', 'off');
    inputPath = fullfile(imagePath,imageFile);
    set(handles.txtInputfilename,'string',inputPath)
    I = imread(inputPath);
    ig = rgb2gray(I);
    
	
    %load image in panel
    axes(handles.axes1);
    hSP = handles.axes1.Parent; %scrollpanel
    hFig = hSP.Parent; %figure
    api = iptgetapi(hSP);
    api.replaceImage(I)
    
    %Panel 2 - show histogram  
    [cell1,cell2]=generateHistogram(I,handles);
    data = struct('img',I, 'inputfile', inputPath, ...
        'ig', ig,'cell1',cell1,'cell2',cell2);
    hObject.UserData = data;
    %clear any previous data
    htable = findobj('Tag','uitableResults');
    set(htable,'Data',[]);
    %or set(htable,'Visible','off')
    %Directions
    updateStatus(handles,'Image loaded: Proceed to set ROI with ROI tool - dbl-click in ROI to save');

% --- Executes on button press in btnAnalysisFiles.
function btnAnalysisFiles_Callback(hObject, eventdata, handles)
% hObject    handle to btnAnalysisFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    fileTypes = {  '*.csv', 'CSV' };
    hCf1 = findobj('Tag', 'editCell1');
    cell1ID=get(hCf1,'String');% 'DS';
    hCf1 = findobj('Tag', 'editCell2');
    cell2ID=get(hCf1,'String');%'SBAC';
    [csvFile, csvPath, filterIndex] = ...
      uigetfile(fileTypes, ...
		'Select analysis files', ...
		'MultiSelect', 'on');
    if iscell(csvFile)
        nbfiles = length(csvFile);
        for i=1:nbfiles
            if (size(strfind(csvFile{i},cell1ID)) > 0)
                cell1 = csvFile{i};
            elseif (size(strfind(csvFile{i}, cell2ID)) > 0)
                cell2 = csvFile{i};
            end
        end
        if (strcmp(cell1, cell1ID) || strcmp(cell2, cell2ID))
            errstr = sprintf('ERROR: Requires 2 CSV files with filenames containing %s and %s',cell1ID,cell2ID);
            updateStatus(handles,errstr);
        else
            analysisfiles = strcat(cell1, ', ', cell2);
            set(handles.editAnalysisfiles,'string',analysisfiles);
            csv1 = fullfile(csvPath, cell1);
            csv2 = fullfile(csvPath, cell2);
            status = strcat('Cell 1: ', csv1,'; Cell 2: ', csv2);
            updateStatus(handles,status);
            %save to userdata
            data = struct('csvPath', csvPath, 'cell1file', cell1, 'cell2file',cell2);
            hObject.UserData = data;
        end
    else
       updateStatus(handles,'ERROR: Requires 2 CSV files'); 
    end

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
    setTitlePlot2(handles,'Histogram');  
    hBg = findobj('Tag','radioWhite');
    bg = 256;
    if (get(hBg,'Value') == 0)
        bg = 0;
    end
    [cell1,cell2] = SplitCells(I,bg);





% --------------------------------------------------------------------
function menu_File_Callback(hObject, eventdata, handles)
% hObject    handle to menu_File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



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
h = findobj('Tag','btnBrowser');
data = h.UserData;
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
    roi_area = polyarea(xi,yi)
    %temp - copy to file
    imwrite(BW, 'tmp_roi.tif');
    %save data - workaround cannot save userdata in radio button??
    dataroi = struct('roi', BW, 'xi',xi,'yi',yi,'x',x,'y',y);
    hObject.UserData = dataroi;
    %show mask image in plot 2
    setTitlePlot2(handles,'');   
    axes(handles.axes2);
    cla;
    imshow(BW)
    title('ROI Mask');
    %Results table
    RCols={'ROI NumberPixels', 'ROIArea'};
    RData=[numPixels, roi_area];
    htable = findobj('Tag','uitableResults');
    set(htable,'data',RData,'ColumnName',RCols);
    %Directions
    status= sprintf('ROI set: Proceed to Identify button')
    updateStatus(handles,status);
else
    hObject.UserData = {};
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


function updateStatus(handles,string)
    if strfind(string,'ERROR')
        set(handles.textOutput, 'string', string, 'Foreground', [1 0 0])
    else
        set(handles.textOutput, 'string', string, 'Foreground', [0 0 0])
    end
    
function setTitlePlot1(handles,string)
    set(handles.txtTitleplot1, 'string', string)
function setTitlePlot2(handles,string)
    set(handles.txtTitleplot2, 'string', string)


% --- Executes on button press in btnIdentify.
function btnIdentify_Callback(hObject, eventdata, handles)
% hObject    handle to btnIdentify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get grayscale img
h = findobj('Tag','btnBrowser');
data = h.UserData;

%show figs
axes(handles.axes2);
%hFig = handles.axes2.Parent;

cla;

%Run analyser
%N = NeuritesAnalyser('example.tif','tmp_roi.tif');
N = NeuritesAnalyser(data.inputfile,'tmp_roi.tif',data.cell1,data.cell2);
%hIm = imshow(data.ig);
%hSP = imscrollpanel(hFig,hIm); % Handle to scroll panel.
%set(hSP,'Units','normalized',...
%        'Position',[0.535 .369 0.421 .582]) %match with GUI box
title('Analysed ROI');
setTitlePlot2(handles,'');
%configuration
h1 = findobj('Tag', 'editMinlength');
minlength = str2double(get(h1,'String'));
h2 = findobj('Tag', 'editTolerance');
tolerance = str2double(get(h2,'String'));
h3 = findobj('Tag', 'checkboxShowplots');
showplots = get(h3,'Value');
types=[];
h4 = findobj('Tag', 'checkboxEnpassant');
if (get(h4,'Value') > 0)
    types(length(types)+1)=1;
end
h5 = findobj('Tag', 'checkboxEndpoint');
if (get(h5,'Value') > 0)
    types(length(types)+1)=2;
end
h6 = findobj('Tag', 'checkboxIntersection');
if (get(h6,'Value') > 0)
    types(length(types)+1)=3;
end

%Run analysis
N = N.findSegments(minlength,tolerance);
N = N.analyseSynapses(showplots,types);
csvdata = handles.btnAnalysisFiles.UserData;
if ( ~isempty(csvdata))
    csv1 = fullfile(csvdata.csvPath, csvdata.cell1file);
    csv2 = fullfile(csvdata.csvPath, csvdata.cell2file);
    N = N.measureSynapses(types,csv1,csv2);
end
ctr = 0;
axes(handles.axes1);
hold on;
mycolors=['m' 'b' 'c'];
mytypes =['En passant' 'End point' 'Intersection'];
RData = [];
%show data in table
colnames = {'Type','Cell 1 X','Cell 1 Y','Cell 2 X','Cell 2 Y', 'Synapse length'};
htable = findobj('Tag','uitableResults');
for i=1:length(N.Synapses)
    syn = N.Synapses{i};
    if (ismember(syn.SynapseType,types))
        ctr = ctr + 1;
        %plot onto original
        plot(syn.MedianC1(:,1), syn.MedianC1(:,2),'color',mycolors(syn.SynapseType),'marker','X','linestyle','none','LineWidth', 2);
        plot(syn.MedianC2(:,1), syn.MedianC2(:,2),'color',mycolors(syn.SynapseType),'marker','X','linestyle','none','LineWidth', 2);
        %show data in table
        row1 = [syn.SynapseType syn.MedianC1 syn.MedianC2 length(syn.RegionC1)];
        RData = cat(1,RData,row1);
    end
end
hold off;
set(htable,'data',RData,'ColumnName',colnames);
%save to file
outputfile = 'tmp_data.csv';
csvwrite(outputfile,RData);
%save analysis
dataroi = struct('N', N, 'numsynapses',ctr);
hObject.UserData = dataroi;
status = sprintf('Found %d synapses. Saved to file:%s', ctr,outputfile);
updateStatus(handles,status);




% --- Executes on button press in btnSplitcells.
function btnSplitcells_Callback(hObject, eventdata, handles)
% hObject    handle to btnSplitcells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Use RGB image
h = findobj('Tag','btnBrowser');
data = h.UserData;
I = data.img;
hBg = findobj('Tag','radioWhite');
bg = 256;
if (get(hBg,'Value') == 0)
    bg = 0;
end
[cell1,cell2] = SplitCells(I,bg);



function editAnalysisfiles_Callback(hObject, eventdata, handles)
% hObject    handle to editAnalysisfiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editAnalysisfiles as text
%        str2double(get(hObject,'String')) returns contents of editAnalysisfiles as a double


% --- Executes during object creation, after setting all properties.
function editAnalysisfiles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editAnalysisfiles (see GCBO)
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


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Results table
numPixels = 15;
roi_area = 30;
colnames = {'ROI NumberPixels','ROIArea'};
RData=[numPixels roi_area];
htable = findobj('Tag','uitableResults');
set(htable,'data',RData,'ColumnName',colnames);



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
