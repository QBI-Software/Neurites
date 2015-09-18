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

% Last Modified by GUIDE v2.5 16-Sep-2015 15:40:33

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
    hIm = imshow('synapse_EM.jpg');
    hSP = imscrollpanel(hFig,hIm); % Handle to scroll panel.
    set(hSP,'Units','pixels',...
        'Position',[50 334 551 451]) %match with GUI box
	    
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
    data = struct('img',I, 'imagePath',imagePath,'inputfile', inputPath, ...
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
    if (iscell(cell1ID))
        cell1ID = strjoin(cell1ID);
    end
    hCf1 = findobj('Tag', 'editCell2');
    cell2ID=get(hCf1,'String');%'SBAC';
    if (iscell(cell2ID))
        cell2ID = strjoin(cell2ID);
    end
    hfiles = findobj('Tag','btnBrowser');
    hfdata = hfiles.UserData;
    csvPath = hfdata.imagePath;
    [csvFile, csvPath, filterIndex] = ...
      uigetfile(fileTypes, ...
		'Select analysis files', ...
		'MultiSelect', 'on',csvPath);
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
            errordlg(errstr,'CSV Files error')
            %updateStatus(handles,errstr);
        else
            analysisfiles = strcat(cell1, ', ', cell2);
            set(handles.editAnalysisfiles,'string',analysisfiles);
            csv1 = fullfile(csvPath, cell1);
            csv2 = fullfile(csvPath, cell2);
            %status = strcat('Cell 1: ', csv1,'; Cell 2: ', csv2);
            status = sprintf('Cell 1: %s\nCell 2: %s', csv1,csv2);
            updateStatus(handles,status);
            %save to userdata
            data = struct('csvPath', csvPath, 'cell1file', cell1, 'cell2file',cell2);
            hObject.UserData = data;
            %load config data if present
            configfile =fullfile(csvPath, 'neurites_config.csv');
            if (exist(configfile, 'file') == 2)
                M = readtable(configfile);
                loadConfig(M,handles);
            end
        end
    else
       errordlg('Requires 2 CSV files','CSV Files error')
    end

function loadConfig(M, handles)
    hCell1 = findobj('Tag','editCell1');
    hCell2 = findobj('Tag','editCell2');
    hScale = findobj('Tag','editScale');
    hSX = findobj('Tag','editShiftx');
    hSY = findobj('Tag','editShifty');
    hFit = findobj('Tag','editFit');
    hMin = findobj('Tag','editMinlength');
    hTol = findobj('Tag','editTolerance');
    set(hCell1,'String',M.Cell1);
    set(hCell2,'String',M.Cell2);
    set(hScale,'String',num2str(M.Scale));
    set(hFit,'String',num2str(M.Fit));
    set(hSX,'String',num2str(M.Shiftx));
    set(hSY,'String',num2str(M.Shifty));
    set(hMin,'String',num2str(M.Minlength));
    set(hTol,'String',num2str(M.Tolerance));

    % --- Executes on button press in pushbutton9.
function saveConfig_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
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
    Cell1 = {get(hCell1,'String')};
    Cell2 = {get(hCell2,'String')};
    Scale = [str2double(get(hScale,'String'))];
    Fit = [str2double(get(hFit,'String'))];
    Shiftx = [str2double(get(hSX,'String'))];
    Shifty = [str2double(get(hSY,'String'))];
    Minlength = [str2double(get(hMin,'String'))];
    Tolerance = str2double(get(hTol,'String'));

    T = table(Cell1, Cell2, Scale, Fit, Shiftx, Shifty, Minlength,Tolerance);
    %save to file
    csvdata = handles.btnAnalysisFiles.UserData;
    outputfile = 'neurites_config.csv';
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
    roi_area = polyarea(xi,yi);
    %temp - copy to file
    hfiles = findobj('Tag','btnBrowser');
    hfdata = hfiles.UserData;
    roifile = fullfile(hfdata.imagePath, 'neurites_roi.tif');
    imwrite(BW, roifile);
    %save data - workaround cannot save userdata in radio button??
    dataroi = struct('roi', BW, 'xi',xi,'yi',yi,'x',x,'y',y);
    hObject.UserData = dataroi;
    %show mask image in plot 2
    %setTitlePlot2(handles,'');   
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


% --- Executes on button press in btnIdentify.
function btnIdentify_Callback(hObject, eventdata, handles)
% hObject    handle to btnIdentify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get grayscale img
h = findobj('Tag','btnBrowser');
data = h.UserData;
hRoi = findobj('Tag','radioROI');
roidata = hRoi.UserData;
csvdata = handles.btnAnalysisFiles.UserData;
if (~isempty(roidata))
    roi = roidata.roi;
    xi = roidata.xi;
    yi = roidata.yi;
elseif(~isempty(csvdata))
    roi = fullfile(csvdata.csvPath, 'neurites_roi.tif');
else
    errordlg('Please set ROI first','ROI error!')
    return
end


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

%configuration
h1 = findobj('Tag', 'editMinlength');
minlength = str2double(get(h1,'String'));
h2 = findobj('Tag', 'editTolerance');
tolerance = str2double(get(h2,'String'));
h3 = findobj('Tag', 'checkboxShowplots');
showplots = get(h3,'Value');
types =[1];
hScale = findobj('Tag','editScale');
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
step = steps;
waitbar(step / steps)
close(hwb);

ctr = 0;
axes(handles.axes1);
hold on;
%Plot ROI for comparison
if (~isempty(roidata))
    plot(xi, yi, '--b','LineWidth', 1);
end
mycolors=['m' 'b' 'c'];
%mytypes =['En passant' 'End point' 'Intersection'];
htable = findobj('Tag','uitableResults');
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
%Save image
%hFig = handles.axes1.Parent;
%hIm = hFig.getimage;
%hSP = imscrollpanel(hFig,hIm)
F=getframe(handles.axes1); %select axes in GUI
figure(); %new figure
image(F.cdata); %show selected axes in new figure
saveas(gcf, fullfile(pathname, 'neurites_img.png')); %save figure
close(gcf); %and close it


%Save data
[colnames,T] = N.generateTable(types,cell1label, cell2label);
set(htable,'data',T,'ColumnName',colnames);
%save to file
outputfile = fullfile(pathname, 'neurites_data.csv');
saveDataFile(outputfile, colnames, T);
%csvwrite(outputfile,T);
%T = cat(1,colnames,T);
%writetable(T,outputfile);
%save analysis
dataN = struct('N', N, 'numsynapses',ctr);
hObject.UserData = dataN;
if (nocsvflag)
    csvmsg ='Load CSV files to get full statistics.';
else
    csvmsg = 'Statistics from CSV files.';
end

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


% --- Executes on button press in btnRegister.
function btnRegister_Callback(hObject, eventdata, handles)
% hObject    handle to btnRegister (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = findobj('Tag','btnBrowser');
data = h.UserData;
I = data.cell1;
csvdata = handles.btnAnalysisFiles.UserData;
hScale = findobj('Tag','editScale');
hSX = findobj('Tag','editShiftx');
hSY = findobj('Tag','editShifty');
if ( ~isempty(csvdata))
    csv1 = fullfile(csvdata.csvPath, csvdata.cell1file);
    T1 = readtable(csv1);
    hscale = str2double(get(hScale,'String'));
    shiftx = str2double(get(hSX,'String'));
    shifty = str2double(get(hSY,'String'));
    %plot
    X1 = (T1.StartX * hscale) + shiftx;
    Y1 = (T1.StartY * hscale) + shifty;
    axes(handles.axes2);
    imshow(I);
    hold on;
    plot(X1,-Y1,'color','c','LineStyle','-','LineWidth', 2);
    hold off;
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


% --- Executes on button press in btnReview.
function btnReview_Callback(hObject, eventdata, handles)
% hObject    handle to btnReview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox('Under development!','Warn');


% --------------------------------------------------------------------
function Menu_About_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_About (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Message = {'QBI Williams Neurites Synapse Analyser' ...
    'Description: Developed for Simon de Croft (Williams lab) for the '...
    'detection and analysis of synaptic regions of two neurons.' ...
    'Developed by Liz Cooper-Williams (e.cooperwilliams@uq.edu.au), (c)Copyright QBI 2015' ...
    'Version: 1.1 (Sep 2015)'};
msgbox(Message,'About')

% --------------------------------------------------------------------
function Menu_Help_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Message = {'Instructions for use' ...
    '1. Load image(tif or jpg), check histogram and split cells' ...
    '2. Load 2 CSV data files for each cell and check names(eg DS.csv and SBAC.csv)' ...
    '3. Configure CSV overlay with Register CSV - adjust configuration values then save config' ...
    '4. Draw ROI with ROI tool and double click in region to save' ...
    '5. Run analysis with Identify synaptic regions' ...
    '6. Check results'};
h = msgbox(Message,'Instructions','Help');


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);
