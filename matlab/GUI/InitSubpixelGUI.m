function handles = InitSubpixelGUI(hObject, handles)

global Image;
Image = [];

handles.currentPanel = 0;
handles.currentFunction = 0;
handles.edges = 0;
handles.orientation = 0;
handles.xCenter = 0;
handles.yCenter = 0;
handles.innerRadius = 0;
handles.outerRadius = 0;
handles.innerIntensity = 0;
handles.outerIntensity = 0;

% set final height of main window
posGUI = get(handles.MainWindow, 'Position');
posGUI(3) = 400;
posGUI(4) = 360;
set(handles.MainWindow, 'Position', posGUI);

% disable other panels
panels = [handles.NoisePanel, handles.DetectionPanel];
set(findall(panels, '-property', 'enable'), 'enable', 'off');
set(handles.checkboxShowEdges, 'enable', 'off');
set(handles.checkboxShowNormals, 'enable', 'off');

% init image figure
iptsetpref('ImshowBorder','tight');
%set(0,'DefaultFigureCloseRequestFcn',@MyCloseFigure);
handles.currentFigure = figure('Visible', 'off', 'Name', 'Image');
end


function MyCloseFigure(src, event)

% access main window
mainWindow = findobj('Tag','MainWindow');
if mainWindow == []
    return;
end

% delete figure
handles = guidata(mainWindow);
delete(src);

% reset images
global Image;
Image = [];

% disable other panels
if (src ~= mainWindow)
    panels = [handles.NoisePanel, handles.DetectionPanel];
    set(findall(panels, '-property', 'enable'), 'enable', 'off');
    if (handles.currentPanel ~= 0)
        set(handles.currentPanel, 'Visible', 'off');
    end
end

% if closing main window, disable the figure
if (src == mainWindow)
    if (ishandle(handles.currentFigure) ~= 0)
        delete(handles.currentFigure);
    end
end
end