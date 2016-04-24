%C1431388 
function varargout = LightGUIWorking(varargin)
% LightGUIWorking MATLAB code for LightGUIWorkingWorking.fig
%      LightGUIWorking, by itself, creates a new LightGUIWorkingWorking or raises the existing
%      singleton*.
%
%      H = LightGUIWorking returns the handle to a new LightGUIWorkingWorking or the handle to
%      the existing singleton*.
%
%      LightGUIWorking('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LightGUIWorking.M with the given input arguments.
%
%      LightGUIWorking('Property','Value',...) creates a new LightGUIWorking or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LightGUIWorking_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LightGUIWorking_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LightGUIWorking

% Last Modified by GUIDE v2.5 17-Dec-2015 17:20:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LightGUIWorking_OpeningFcn, ...
                   'gui_OutputFcn',  @LightGUIWorking_OutputFcn, ...
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


% --- Executes just before LightGUIWorking is made visible.
function LightGUIWorking_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LightGUIWorking (see VARARGIN)

% Choose default command line output for LightGUIWorking
handles.output = hObject;

%instigate all gloabal variables that will be needed
global mainlight mirrorlist mirrorcounter lightLines reflections reflectioncounter linecounter lines intensityCounter

%create 5 blank structure arrays, mainlight holds original light beam
%plots, mirrorlist holds all placed mirror plots, reflections holds all the
%reflection plots, lines and lightlines holds the line objects rather than
%the x,y plots of the lightbeams.
mainlight = struct('x1', {}, 'x2', {}, 'y1', {}, 'y2', {});
mirrorlist = struct('x1', {}, 'x2', {}, 'y1', {}, 'y2', {}, 'isCirlce', {});
reflections = struct('x1', {}, 'x2', {}, 'y1', {}, 'y2', {}, 'Intensity', {});
lightLines = struct('lineobject', {});
lines = struct('line', {});

%4 counters that are used to keep track of the number of mirrors placed,
%reflections that occur and the intensity of the light.
mirrorcounter = 0;
reflectioncounter = 0;
linecounter = 0;
intensityCounter = 100;

%preparing the axes to be plotted on. by plotting 4 points in the corners
%of the axis it stop it scaling when an object is placed (along with the
%scaling properties in GUIDE.
axes(handles.PlotArea)
guidata(hObject, handles);
hold on
plot(0,0)
plot(0,1)
plot(1,0)
plot(1,1)

%set all buttons except load and plot light beam disabled.
set(handles.MirrorStraight, 'Enable', 'off')
set(handles.MirrorCurve, 'Enable', 'off')
set(handles.Save, 'Enable', 'off')
set(handles.Reset, 'Enable', 'off')
set(handles.Delete, 'Enable', 'off')

% UIWAIT makes LightGUIWorking wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LightGUIWorking_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function updateLight(thisLight, isMain)
%updates light beam when new mirror is placed
%instigate gloabal variables needed for this function.
global mirrorlist mirrorcounter reflections reflectioncounter linecounter lines intensityCounter
    %create variables needed to calculate the closest mirror placed 
    bestMirror = 0;
    bestDistance = inf;
    bestsectx = 0;
    bestsecty = 0;
    
    %check if update was called from placing a lightbeam or mirror
    if isMain == 1
        %clear the axes
        cla('reset')
        grid on
        hold on
        axis manual
        %reset reflection counter
        reflectioncounter = 0;
        %for each mirror in mirrorlist
        for i = 1: 1: mirrorcounter,
            %if it is a straight mirror
            if mirrorlist(i).isCircle == 0
                %plot each mirror using the 2 x,y coordinates and the line
                %function
                mirrorx = [mirrorlist(i).x1, mirrorlist(i).x2];
                mirrory = [mirrorlist(i).y1, mirrorlist(i).y2];
                line(mirrorx, mirrory);
                %else if it is a circular mirror
            else
                %plot circluar mirror using the center point, the radius
                %and 360 degree angle.
                centerpoint = [mirrorlist(i).x1, mirrorlist(i).y1];
                radius = mirrorlist(i).x2;
                angle = 1:0.001:360;
                xPoint = (radius*cos(angle) + centerpoint(1));
                yPoint = (radius*sin(angle) + centerpoint(2));
                plot(xPoint, yPoint)
            end
        end
    end
    %calculates the vector of a longer light beam than plotted, then displaying that
    %light, saving light beam co ordinates in mainlight. this is so the
    %light beam goes to the edge of the axes
    lightvec = [(thisLight(1).x2 - thisLight(1).x1), (thisLight(1).y2 - thisLight(1).y1)];
    lightdistance = sqrt(lightvec(1)*lightvec(1) + lightvec(2)*lightvec(2));
    normalLight = [(2*(lightvec(1)/lightdistance)),(2*(lightvec(2)/lightdistance))];
    theLightx = [(thisLight(1).x1), (thisLight(1).x1 + normalLight(1))];
    theLighty = [(thisLight(1).y1), (thisLight(1).y1 + normalLight(2))];
    thisLight(1).x2 = theLightx(2);
    thisLight(1).y2 = theLighty(2);
    %for each mirror plotted
    for i = 1: 1: mirrorcounter;
        %check if its a straight mirror
        if mirrorlist(i).isCircle == 0
            %find the intercepts of the mirror and the light
            mirrorx = [mirrorlist(i).x1, mirrorlist(i).x2];
            mirrory = [mirrorlist(i).y1, mirrorlist(i).y2];
            [sectx, secty] = polyxpoly([thisLight.x1, thisLight.x2], [thisLight.y1, thisLight.y2], mirrorx, mirrory);
            %if there is an intersect found
            if isempty(sectx) == 0
                %for each intersect (if many) find the closest one to the
                %light source and save that as best intersect
                for a = 1: 1: length(sectx),
                    lightvectorx = (sectx(a) - thisLight(1).x1);
                    lightvectory = (secty(a) - thisLight(1).y1);
                    distance = sqrt((lightvectorx*lightvectorx) + (lightvectory*lightvectory));
                    if distance < bestDistance && distance > 0.00000001
                        bestDistance = distance;
                        bestMirror = i;
                        bestsectx = sectx;
                        bestsecty = secty;
                    end
                end
             end
        else
            %else for circular mirrors find the intersects and store best
            %intersect
            centerpoint = [mirrorlist(i).x1, mirrorlist(i).y1];
            radius = mirrorlist(i).x2;
            for angle = 1:0.1:360, 
                xPoint = (radius*cos(angle) + centerpoint(1));
                yPoint = (radius*sin(angle) + centerpoint(2));
                sectlinex = [xPoint, (xPoint+0.001)];
                sectliney = [yPoint, (yPoint+0.001)];
                [sectx, secty] = polyxpoly([thisLight.x1, thisLight.x2], [thisLight.y1, thisLight.y2], sectlinex, sectliney); 
                if isempty(sectx) == 0
                    for a = 1: 1: length(sectx),
                        lightvectorx = (sectx(a) - thisLight(1).x1);
                        lightvectory = (secty(a) - thisLight(1).y1);
                        distance = sqrt((lightvectorx*lightvectorx) + (lightvectory*lightvectory));
                        if distance < bestDistance && distance > 0.00000001
                            bestDistance = distance;
                            bestMirror = i;
                            bestsectx = sectx;
                            bestsecty = secty;
                        end
                    end
                end
            end
        end
    end
    %assign best x ,y intersect to sectx, secty   
    sectx = bestsectx;
    secty = bestsecty;
    %if a mirror is found to intersect the light beam
    if bestMirror ~= 0
        %adjust the end of the light beam to stop at the intersect
        thisLight(1).x2 = sectx;
        thisLight(1).y2 = secty;
        normalvec = [];
        %if the intersect mirror is not a circle
        if mirrorlist(bestMirror).isCircle == 0
            %calculate the unit vector of the  mirror
            mirrorx = mirrorlist(bestMirror).x1 - mirrorlist(bestMirror).x2;
            mirrory = mirrorlist(bestMirror).y1 - mirrorlist(bestMirror).y2;
            mirrorvec = [mirrorx, mirrory];
            %calculate the normal of the mirror
            normalvec = [mirrorvec(2), (-1*mirrorvec(1))];
            %calculate the unit vector of the light beam
            lightvec = [(sectx - thisLight(1).x1), (secty - thisLight(1).y1)];
            %check to see which side of the mirror the light is coming from
            %and adjust the normal accordingly
            if dot(normalvec, lightvec) < 0 
                normalvec = [(-1*mirrorvec(2)),mirrorvec(1)]; 
            end
        else
            %else if the mirror is a circle, note the centre point,
            %calculate the normal and normalise the normal line
            centerPoint = [mirrorlist(bestMirror).x1, mirrorlist(bestMirror).y1];
            normalvec = [(sectx - mirrorlist(bestMirror).x1), (secty - mirrorlist(bestMirror).y1)];
            distance = sqrt(normalvec(1)*normalvec(1) + normalvec(2)*normalvec(2));
            normalvec = [(normalvec(1)/distance), (normalvec(2)/distance)];
        end
        
        distance = sqrt((normalvec(1)*normalvec(1)) + (normalvec(2)*normalvec(2)));
        normnormal = [(normalvec(1)/distance), (normalvec(2)/distance)];
        %calculate the reflective line 
        refx = (2*dot(lightvec, normnormal)*normnormal(1));
        refy = (2*dot(lightvec, normnormal)*normnormal(2));
        rx = lightvec(1) - refx;
        ry = lightvec(2) - refy;
        %count reflection, add reflection co ordinates to reflections
        %structure
        reflectioncounter = reflectioncounter + 1;
        reflections(reflectioncounter).x1 = sectx;
        reflections(reflectioncounter).y1 = secty;
        reflections(reflectioncounter).x2 = sectx + rx;
        reflections(reflectioncounter).y2 = secty + ry;
        reflections(reflectioncounter).Intensity = intensityCounter;
        %check intensity of the light, if it is too low, delete last mirror
        %placed and display message box warning the user.
        if reflectioncounter == 10
            Delete_Callback()
            msgbox('Light Intensity is too Low, Cannot Reflect Light.', 'Intensity', 'warn')
        end
        %recursive call, instead of passing the main light, you pass the
        %reflection line as the main light.
        updateLight(reflections(reflectioncounter), 0)
        
    end
    %once all calculations are complete plot the lght line and all its
    %reflections
    linecounter = linecounter + 1;
    lines(linecounter).line = line([thisLight(1).x1, thisLight(1).x2], [thisLight(1).y1, thisLight(1).y2], 'Color', 'y');
    
    
    

% --- Executes on button press in LightBeam.
% This function allows the user to plot the inital light beam
function LightBeam_Callback(hObject, eventdata, handles)
% hObject    handle to LightBeam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)  
%intigate the gloabal variables needed for the function
global mainlight
axes(handles.PlotArea)
%get 2 input points from the user as assign the x and y values to the
%variables lightx and lighty
[lightx, lighty] = ginput(2);
%assign the x and y values to the main light structure
mainlight(1).x1 = lightx(1);
mainlight(1).x2 = lightx(2);
mainlight(1).y1 = lighty(1);
mainlight(1).y2 = lighty(2);
%call the recursive function update light, passing it the mainlight
%structure
updateLight(mainlight, 1)
%disable the lightbeam and load buttons, and enable all the rest
set(handles.LightBeam, 'Enable', 'off')
set(handles.Load, 'Enable', 'off')
set(handles.MirrorStraight, 'Enable', 'on')
set(handles.MirrorCurve, 'Enable', 'on')
set(handles.Save, 'Enable', 'on')
set(handles.Reset, 'Enable', 'on')



% --- Executes on button press in MirrorStraight.
%this function allows the user to plot a straigh mirror
function MirrorStraight_Callback(hObject, eventdata, handles)
% hObject    handle to MirrorStraight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%instigate gloabl variables need for the function
global mainlight mirrorlist mirrorcounter
    %get 2 input points from the user as assign the x and y values to the
    %variables smirrorx and smirrory
    [smirrorx, smirrory] = ginput(2);
    %increment the mirror counter
    mirrorcounter = mirrorcounter + 1;
    %assign mirror coordinates the the mirrorlist structure with the index
    %of mirror counter
    mirrorlist(mirrorcounter).x1 = smirrorx(1);
    mirrorlist(mirrorcounter).x2 = smirrorx(2);
    mirrorlist(mirrorcounter).y1 = smirrory(1);
    mirrorlist(mirrorcounter).y2 = smirrory(2);
    %it is not a circular mirror
    mirrorlist(mirrorcounter).isCircle = 0;
    %call recurisve function passing it the mian light
    updateLight(mainlight, 1);
    %enable delete function, as there is mirrors that can be deleted now. 
    set(handles.Delete, 'Enable', 'on')

% --- Executes on button press in MirrorCurve.
%allows the user to plot a circular mirror. 
function MirrorCurve_Callback(hObject, eventdata, handles)
% hObject    handle to MirrorCurve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%instigate gloabl variables needed for the function
global mirrorcounter mainlight mirrorlist
    %get 2 input points from the user as assign the x and y values to the
    %variables circlex and circley where first point is the centre point
    %and the second click is the radius
    [circlex, circley] = ginput(2);
    %calculate the radius of the circle
    circlevec = [(circlex(2) - circlex(1)), (circley(2) - circley(1))];
    radius = sqrt(circlevec(1)*circlevec(1) + circlevec(2)*circlevec(2));
    %increment the mirror counter
    mirrorcounter = mirrorcounter + 1;
    %assign the values of the circle to the mirrorlist structure with the
    %index of mirror counter.
    mirrorlist(mirrorcounter).x1 = circlex(1);
    mirrorlist(mirrorcounter).x2 = radius;
    mirrorlist(mirrorcounter).y1 = circley(1);
    mirrorlist(mirrorcounter).y2 = 0;
    %it is a circle
    mirrorlist(mirrorcounter).isCircle = 1;
    %call recursive function passing it main light
    updateLight(mainlight, 1);
    %enable the delete function.
    set(handles.Delete, 'Enable', 'on')
    
% --- Executes on button press in Load.
%allows the user to load a matlab figure.
function Load_Callback(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%prompt the user to type a file name to be loaded
name = inputdlg('Input File Name:', 'Load Figure');
[file] = name{1:1};
%open file.fig
openfig(file)

% --- Executes on button press in Save.
%allows the user to save the current experiment
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%promts the user to input a name to save the figure as
name = inputdlg('Input File Name:', 'Save Figure');
[file] = name{1:1};
%save file.fig
savefig(file)

% --- Executes on button press in Reset.
%resets the axes and functions
function Reset_Callback(hObject, eventdata, handles)
% hObject    handle to Reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%instigates all globl variables needed for the function
global lines mirrorcounter linecounter reflectioncounter reflections mainlight mirrorlist lightLines
%clear and reset that axes
cla('reset')
grid on
hold on
axis manual
%reset/override all of the structures/variables created at the start of the
%program
mainlight = struct('x1', {}, 'x2', {}, 'y1', {}, 'y2', {});
mirrorlist = struct('x1', {}, 'x2', {}, 'y1', {}, 'y2', {}, 'isCirlce', {});
reflections = struct('x1', {}, 'x2', {}, 'y1', {}, 'y2', {});
lightLines = struct('lineobject', {});
lines = struct('line', {});

mirrorcounter = 0;
reflectioncounter = 0;
linecounter = 0;

guidata(hObject, handles);
hold on
%replot the corner points to stop scaling.
plot(0,0)
plot(0,1)
plot(1,0)
plot(1,1)
%reset all the buttons
set(handles.LightBeam, 'Enable', 'on')
set(handles.MirrorStraight, 'Enable', 'off')
set(handles.MirrorCurve, 'Enable', 'off')
set(handles.Save, 'Enable', 'off')
set(handles.Load, 'Enable', 'on')
set(handles.Reset, 'Enable', 'off')
set(handles.Delete, 'Enable', 'off')

% --- Executes on button press in Delete.
function Delete_Callback(hObject, eventdata, handles)
% hObject    handle to Delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global mirrorcounter mainlight
mirrorcounter = mirrorcounter - 1;
updateLight(mainlight, 1);
