function varargout = neurites_review(varargin)
% NEURITES_REVIEW MATLAB code for neurites_review.fig
%      NEURITES_REVIEW, by itself, creates a new NEURITES_REVIEW or raises the existing
%      singleton*.
%
%      H = NEURITES_REVIEW returns the handle to a new NEURITES_REVIEW or the handle to
%      the existing singleton*.
%
%      NEURITES_REVIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NEURITES_REVIEW.M with the given input arguments.
%
%      NEURITES_REVIEW('Property','Value',...) creates a new NEURITES_REVIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before neurites_review_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to neurites_review_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help neurites_review

% Last Modified by GUIDE v2.5 16-Sep-2015 16:18:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @neurites_review_OpeningFcn, ...
                   'gui_OutputFcn',  @neurites_review_OutputFcn, ...
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


% --- Executes just before neurites_review is made visible.
function neurites_review_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to neurites_review (see VARARGIN)

% Choose default command line output for neurites_review
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes neurites_review wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = neurites_review_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btnBack.
function btnBack_Callback(hObject, eventdata, handles)
% hObject    handle to btnBack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkboxAccept.
function checkboxAccept_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxAccept (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxAccept


% --- Executes on button press in btnNext.
function btnNext_Callback(hObject, eventdata, handles)
% hObject    handle to btnNext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function editXY_Callback(hObject, eventdata, handles)
% hObject    handle to editXY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editXY as text
%        str2double(get(hObject,'String')) returns contents of editXY as a double


% --- Executes during object creation, after setting all properties.
function editXY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editXY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editC1label_Callback(hObject, eventdata, handles)
% hObject    handle to editC1label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editC1label as text
%        str2double(get(hObject,'String')) returns contents of editC1label as a double


% --- Executes during object creation, after setting all properties.
function editC1label_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editC1label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editSomaC1_Callback(hObject, eventdata, handles)
% hObject    handle to editSomaC1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSomaC1 as text
%        str2double(get(hObject,'String')) returns contents of editSomaC1 as a double


% --- Executes during object creation, after setting all properties.
function editSomaC1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSomaC1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editEndC1_Callback(hObject, eventdata, handles)
% hObject    handle to editEndC1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editEndC1 as text
%        str2double(get(hObject,'String')) returns contents of editEndC1 as a double


% --- Executes during object creation, after setting all properties.
function editEndC1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editEndC1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editC2label_Callback(hObject, eventdata, handles)
% hObject    handle to editC2label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editC2label as text
%        str2double(get(hObject,'String')) returns contents of editC2label as a double


% --- Executes during object creation, after setting all properties.
function editC2label_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editC2label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editSomaC2_Callback(hObject, eventdata, handles)
% hObject    handle to editSomaC2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSomaC2 as text
%        str2double(get(hObject,'String')) returns contents of editSomaC2 as a double


% --- Executes during object creation, after setting all properties.
function editSomaC2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSomaC2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editEndC2_Callback(hObject, eventdata, handles)
% hObject    handle to editEndC2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editEndC2 as text
%        str2double(get(hObject,'String')) returns contents of editEndC2 as a double


% --- Executes during object creation, after setting all properties.
function editEndC2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editEndC2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
