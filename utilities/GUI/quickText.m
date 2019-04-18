function varargout = quickText(varargin)
% QUICKTEXT MATLAB code for quickText.fig
%      QUICKTEXT, by itself, creates a new QUICKTEXT or raises the existing
%      singleton*.
%
%      H = QUICKTEXT returns the handle to a new QUICKTEXT or the handle to
%      the existing singleton*.
%
%      QUICKTEXT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in QUICKTEXT.M with the given input arguments.
%
%      QUICKTEXT('Property','Value',...) creates a new QUICKTEXT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before quickText_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to quickText_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help quickText

% Last Modified by GUIDE v2.5 08-Oct-2017 17:15:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @quickText_OpeningFcn, ...
                   'gui_OutputFcn',  @quickText_OutputFcn, ...
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


% --- Executes just before quickText is made visible.
function quickText_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to quickText (see VARARGIN)

% Choose default command line output for quickText
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes quickText wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = quickText_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function MainTextArea_Callback(hObject, eventdata, handles)
% hObject    handle to MainTextArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MainTextArea as text
%        str2double(get(hObject,'String')) returns contents of MainTextArea as a double


% --- Executes during object creation, after setting all properties.
function MainTextArea_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MainTextArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on MainTextArea and none of its controls.
function MainTextArea_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to MainTextArea (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if(strcmp(eventdata.Key,'return'))
    mx=get(hObject, 'Max');
    str=get(hObject,'String');
    cols=size(str,2);
    set(hObject,'Max', mx+1);
    set(hObject,'String', [str sprintf('\n')]);
end
