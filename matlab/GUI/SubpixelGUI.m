function varargout = SubpixelGUI(varargin)
% SUBPIXELGUI M-file for SubpixelGUI.fig
%      SUBPIXELGUI, by itself, creates a new SUBPIXELGUI or raises the existing
%      singleton*.
%
%      H = SUBPIXELGUI returns the handle to a new SUBPIXELGUI or the handle to
%      the existing singleton*.
%
%      SUBPIXELGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SUBPIXELGUI.M with the given input arguments.
%
%      SUBPIXELGUI('Property','Value',...) creates a new SUBPIXELGUI or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SubpixelGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SubpixelGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SubpixelGUI

% Last Modified by GUIDE v2.5 14-Oct-2015 14:19:52


%******************************************************
%// FOLLOWIN CODE IS GENERATED AUTOMATIALLY BY GUIDE.
%// SUBPIXEL GUI CODE STARTS NEAR THE END OF THE FILE
%*****************************************************

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SubpixelGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SubpixelGUI_OutputFcn, ...
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

% --- Executes just before SubpixelGUI is made visible.
function SubpixelGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SubpixelGUI (see VARARGIN)

% Choose default command line output for SubpixelGUI
handles.output = hObject;

% call my init function
handles = InitSubpixelGUI(hObject, handles);
guidata(hObject, handles);

% standard initialization
initialize_gui(hObject, handles, false);

% UIWAIT makes SubpixelGUI wait for user response (see UIRESUME)
% uiwait(handles.MainWindow);


% --- Outputs from this function are returned to the command line.
function varargout = SubpixelGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function density_CreateFcn(hObject, eventdata, handles)
% hObject    handle to density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function density_Callback(hObject, ~, handles)
% hObject    handle to density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of density as text
%        str2double(get(hObject,'String')) returns contents of density as a double
density = str2double(get(hObject, 'String'));
if isnan(density)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new density value
handles.metricdata.density = density;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function volume_CreateFcn(hObject, eventdata, handles)
% hObject    handle to volume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function volume_Callback(hObject, eventdata, handles)
% hObject    handle to volume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of volume as text
%        str2double(get(hObject,'String')) returns contents of volume as a double
volume = str2double(get(hObject, 'String'));
if isnan(volume)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new volume value
handles.metricdata.volume = volume;
guidata(hObject,handles);

% --- Executes on button press in calculate.
function calculate_Callback(hObject, eventdata, handles)
% hObject    handle to calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%mass = handles.metricdata.density * handles.metricdata.volume;
%set(handles.mass, 'String', mass);
fprintf('leyendo valores del gui\n');


% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

initialize_gui(gcbf, handles, true);

% --- Executes when selected object changed in CirclePanel.
function CirclePanel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in CirclePanel 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (hObject == handles.english)
    set(handles.text4, 'String', 'lb/cu.in');
    set(handles.text5, 'String', 'cu.in');
    set(handles.text6, 'String', 'lb');
else
    set(handles.text4, 'String', 'kg/cu.m');
    set(handles.text5, 'String', 'cu.m');
    set(handles.text6, 'String', 'kg');
end

% --------------------------------------------------------------------
function initialize_gui(fig_handle, handles, isreset)
% If the metricdata field is present and the reset flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to reset the data.
if isfield(handles, 'metricdata') && ~isreset
    return;
end

handles.metricdata.density = 0;
handles.metricdata.volume  = 0;

%set(handles.density, 'String', handles.metricdata.density);
%set(handles.volume,  'String', handles.metricdata.volume);
%set(handles.mass, 'String', 0);

%set(handles.text4, 'String', 'lb/cu.in');
%set(handles.text5, 'String', 'cu.in');
%set(handles.text6, 'String', 'lb');

% Update handles structure
guidata(handles.MainWindow, handles);


% --- Executes when MainWindow is resized.
function MainWindow_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to MainWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when CirclePanel is resized.
function CirclePanel_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to CirclePanel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over CircleButton.
function CircleButton_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to CircleButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over calculate.
function calculate_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit24 as text
%        str2double(get(hObject,'String')) returns contents of edit24 as a double


% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%*****************************************
%// HERE START SUBPIXEL GUI CODE
%****************************************

function DisplayImage(image, name, resetAspect, handles)
if (ishandle(handles.currentFigure) == 0)
    figure(handles.currentFigure);
end
set(handles.currentFigure, 'Visible', 'on', 'Name', name);
if (resetAspect)
    imshow(image/255);
else
    limits = get(gca, {'xlim','ylim'});  
    imshow(image/255, 'InitialMagnification', 'fit');
    zoom reset;
    set(gca, {'xlim','ylim'}, limits)
end


function createRamp(hObject, handles)

% reading parameters
xCenter = 200 + str2num(get(handles.edit20, 'String'));
yCenter = 200 + str2num(get(handles.edit21, 'String'));
size = 400;
orientation = str2num(get(handles.edit24, 'String')); 
innerIntensity = str2num(get(handles.edit22, 'String'));
outerIntensity = str2num(get(handles.edit23, 'String'));
if (innerIntensity < 0)
    innerIntensity = 0;
    set(handles.edit22, 'String', '0.0');
end
if (innerIntensity > 255)
    innerIntensity = 255;
    set(handles.edit22, 'String', '255.0');
end
if (outerIntensity < 0)
    outerIntensity = 0;
    set(handles.edit23, 'String', '0.0');
end
if (outerIntensity > 255)
    outerIntensity = 255;
    set(handles.edit23, 'String', '255.0');
end

% generate ramp
global Image;
Image = ramp(size, size, xCenter, yCenter, ...
    orientation, innerIntensity, outerIntensity);

% save parameters in handles
global edges;
edges = EdgePixel;
handles.orientation = orientation;
handles.xCenter = xCenter;
handles.yCenter = yCenter;
handles.innerRadius = -9999;  % means infinite
handles.outerRadius = -9999;
handles.innerIntensity = innerIntensity;
handles.outerIntensity = outerIntensity;
guidata(hObject, handles);

% display in figure
DisplayImage(Image, 'Synthetic Ramp', true, handles);

% enable and disable other panels
panels = [handles.NoisePanel, handles.DetectionPanel];
set(findall(panels, '-property', 'enable'), 'enable', 'on');
set(handles.checkboxShowEdges, 'enable', 'off');
set(handles.checkboxShowNormals, 'enable', 'off');


function createCircle(hObject, handles)

% reading parameters
innerRadius = str2num(get(handles.edit19, 'String'));
outerRadius = str2num(get(handles.edit9, 'String'));
margin = 10;
center = margin + outerRadius + 1;
if (center < 200)
    center = 200;
end
size = 2*center - 1;
xCenter = center + str2num(get(handles.edit25, 'String'));
yCenter = center + str2num(get(handles.edit26, 'String'));
innerIntensity = str2num(get(handles.edit10, 'String'));
outerIntensity = str2num(get(handles.edit11, 'String'));
if (innerIntensity < 0)
    innerIntensity = 0;
    set(handles.edit10, 'String', '0.0');
end
if (innerIntensity > 255)
    innerIntensity = 255;
    set(handles.edit10, 'String', '255.0');
end
if (outerIntensity < 0)
    outerIntensity = 0;
    set(handles.edit11, 'String', '0.0');
end
if (outerIntensity > 255)
    outerIntensity = 255;
    set(handles.edit11, 'String', '255.0');
end

% generate circle or ring
global Image;
if (innerRadius < 0.01)
    Image = circle(size, size, xCenter, yCenter, outerRadius, ...
        innerIntensity, outerIntensity, 100);
else
    Image = ring(size, size, xCenter, yCenter, innerRadius, outerRadius, ...
        innerIntensity, outerIntensity, 100);   
end

% save parameters in handles
global edges;
edges = EdgePixel;
handles.xCenter = xCenter;
handles.yCenter = yCenter;
handles.innerRadius = innerRadius;
handles.outerRadius = outerRadius;
handles.innerIntensity = innerIntensity;
handles.outerIntensity = outerIntensity;
guidata(hObject, handles);

% display in figure
DisplayImage(Image, 'Synthetic Circle', true, handles);

% enable and disable other panels
panels = [handles.NoisePanel, handles.DetectionPanel];
set(findall(panels, '-property', 'enable'), 'enable', 'on');
set(handles.checkboxShowEdges, 'enable', 'off');
set(handles.checkboxShowNormals, 'enable', 'off');


function addNoise(hObject, handles)

% reading parameters
percent = str2num(get(handles.edit12, 'String'));
if (percent < 0)
    percent = 0;
    set(handles.edit12, 'String', '0.0');
end

% add noise
global Image;
Image = noise(Image, percent);

% reset edges
global edges;
edges = EdgePixel;
guidata(hObject, handles);

% display in figure
DisplayImage(Image, 'Noisy Image', false, handles);

% enable and disable other panels
set(handles.checkboxShowEdges, 'enable', 'off');
set(handles.checkboxShowNormals, 'enable', 'off');


function subpixelDetector(hObject, handles)
threshold = str2num(get(handles.edit13, 'String'));
if (threshold < 0)
    threshold = 0.01;
    set(handles.edit13, 'String', '0.01');
end
order = str2num(get(handles.edit14, 'String'));
if (order <1.5)
    order = 1;
    set(handles.edit14, 'String', '1');
else
    order = 2;
    set(handles.edit14, 'String', '2');
end
iter = floor(str2num(get(handles.edit18, 'String')));

% detect edges 
global Image;
global edges;
global method;
global restoredImage;
method = get(handles.listbox1, 'Value');
switch method
    case 1
        fprintf('Subpixel detection - basic detector ');
    case 2
        fprintf('Subpixel detection - smooth detector ');
    case 3
        fprintf('Subpixel detection - nearby edge detector ');
    case 4
        fprintf('Subpixel detection - final detector iter 0 ');
    case 5
        fprintf('Subpixel detection - final detector iter 1 ');
    case 6
        fprintf('Subpixel detection - final detector iter n ');
end
fprintf('(threshold=%.1f, order=%d) ...\n', threshold, order);

% call detector
tic;
switch method
    case 1
        edges = basicDetector(Image, threshold, order);
    case 2
        edges = smoothDetector(Image, threshold, order);
    case 3
        edges = nearbyEdgeDetector(Image, threshold, order);
    case 4
        edges = finalDetectorIter0(Image, threshold, order);
    case 5
        edges = finalDetectorIter1(Image, threshold, order);
    case 6 
        restoredImage = Image;
        for n=1:iter
            fprintf ('Iteration %d / %d ...\n', n, iter);
            [edges, restoredImage] = finalDetectorIterN(restoredImage, ...
                threshold, order);
        end
end
toc;

n = size(edges.position, 1);
if (n == 0)
    fprintf('No edge pixels detected. Try to decrease threshold\n');
else
    fprintf('Number of edge pixels detected = %d\n', n);
end

% display edges
showRestoredImage = get(handles.checkboxShowRestoredImage, 'Value');
if method==6 && showRestoredImage
    DisplayImage(restoredImage, 'Edges detected', false, handles);
else
    DisplayImage(Image, 'Edges detected', false, handles);
end
hold on;
seg = 0.6;
if (get(handles.checkboxShowEdges, 'Value') == 1)
    quiver(edges.x-seg/2*edges.ny, edges.y+seg/2*edges.nx, ...
        seg*edges.ny, -seg*edges.nx, 0, 'r.');
end
if (get(handles.checkboxShowNormals, 'Value') == 1)
    quiver(edges.x, edges.y, edges.nx, edges.ny, 0, 'b');
end
hold off;

% in case of synthetic image, show statistics
if (handles.innerRadius < 0)
    statRamp(handles.xCenter, handles.yCenter, handles.orientation, ...
        handles.innerIntensity, handles.outerIntensity, edges);
elseif (handles.innerRadius < 0.01)
    statCircle(handles.xCenter, handles.yCenter, handles.outerRadius, ...
        handles.innerIntensity, handles.outerIntensity, edges);
else
    statRing(handles.xCenter, handles.yCenter, ...
        handles.innerRadius, handles.outerRadius, ...
        handles.innerIntensity, handles.outerIntensity, edges);
end

% enable and disable other panels
set(handles.checkboxShowEdges, 'enable', 'on');
set(handles.checkboxShowNormals, 'enable', 'on');
if method == 6
    set(handles.checkboxShowRestoredImage, 'enable', 'on');
else
    set(handles.checkboxShowRestoredImage, 'enable', 'off');
end

% --- Executes on button press in all the function panels.
function PanelButton_Callback(hObject, eventdata, handles)
% hObject    handle to CircleButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.currentFunction(hObject, handles);


function activePanel(hObject, handles, panel, func)
if (handles.currentPanel ~= 0)
    set(handles.currentPanel, 'Visible', 'off');
end
handles.currentPanel = panel;
handles.currentFunction = func;
guidata(hObject, handles);
pos = get(panel, 'Position');
pos(1) = 190;
posGUI = get(handles.MainWindow, 'Position');
vertMargin = 10;
pos(2) = posGUI(4) - pos(4) - vertMargin;
set(panel, 'Position', pos);
set(panel, 'Visible', 'on');


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
activePanel(hObject, handles, handles.CirclePanel, @createCircle);


% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
activePanel(hObject, handles, handles.RampPanel, @createRamp);


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
activePanel(hObject, handles, handles.GaussianNoisePanel, @addNoise);


% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
activePanel(hObject, handles, handles.DetectorPanel, @subpixelDetector);


% --- Executes on key press 
function EditPanel_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if (strcmp(eventdata.Key, 'return') && handles.currentFigure ~= 0)
    figure(handles.currentFigure);  % to update edit text
    handles.currentFunction(hObject, handles);
    uicontrol(hObject);
end


% --- Executes on button press in checkboxShowEdges.
function checkboxPanel_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxShowEdges (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.currentFigure ~= 0)
    figure(handles.currentFigure);  % to update edit text
    handles.currentFunction(hObject, handles);
end

function updateFigure_Callback(hObject, eventdata, handles)
% this callback is called for all the checkboxes in Detector Panel
% hObject    handle to checkboxShowEdges (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Image;
global edges;
global restoredImage;
global method;
showRestoredImage = get(handles.checkboxShowRestoredImage, 'Value');
if method==6 && showRestoredImage
    DisplayImage(restoredImage, 'Edges detected and restored image', ...
        false, handles);
else
    DisplayImage(Image, 'Edges detected', false, handles);
end
hold on;
seg = 0.6;
if (get(handles.checkboxShowEdges, 'Value') == 1)
    quiver(edges.x-seg/2*edges.ny, edges.y+seg/2*edges.nx, ...
        seg*edges.ny, -seg*edges.nx, 0, 'r.');
end
if (get(handles.checkboxShowNormals, 'Value') == 1)
    quiver(edges.x, edges.y, edges.nx, edges.ny, 0, 'b');
end
hold off;

function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double


% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit26_Callback(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit26 as text
%        str2double(get(hObject,'String')) returns contents of edit26 as a double


% --- Executes during object creation, after setting all properties.
function edit26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
