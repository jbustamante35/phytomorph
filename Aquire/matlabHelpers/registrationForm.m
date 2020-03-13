function varargout = registrationForm(varargin)
% REGISTRATIONFORM MATLAB code for registrationForm.fig
%      REGISTRATIONFORM, by itself, creates a new REGISTRATIONFORM or raises the existing
%      singleton*.
%
%      H = REGISTRATIONFORM returns the handle to a new REGISTRATIONFORM or the handle to
%      the existing singleton*.
%
%      REGISTRATIONFORM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REGISTRATIONFORM.M with the given input arguments.
%
%      REGISTRATIONFORM('Property','Value',...) creates a new REGISTRATIONFORM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before registrationForm_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to registrationForm_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help registrationForm

% Last Modified by GUIDE v2.5 31-Jan-2020 11:16:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @registrationForm_OpeningFcn, ...
                   'gui_OutputFcn',  @registrationForm_OutputFcn, ...
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


% --- Executes just before registrationForm is made visible.
function registrationForm_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to registrationForm (see VARARGIN)

% Choose default command line output for registrationForm
handles.output = hObject;

handles.cyversePW = [];

% add the logo for phytomorph
logo = double(imread('~/phytoMorphTK/logo.png'))/255;
axes(handles.axes1);
imshow(logo,[]);

% add the cyverse logo
logo2 = double(imread('~/phytoMorphTK/cyverse.png'))/255;
axes(handles.axes2);
imshow(logo2,[]);
drawnow

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes registrationForm wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = registrationForm_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function name_Callback(hObject, eventdata, handles)
% hObject    handle to name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of name as text
%        str2double(get(hObject,'String')) returns contents of name as a double


% --- Executes during object creation, after setting all properties.
function name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function email_Callback(hObject, eventdata, handles)
% hObject    handle to email (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of email as text
%        str2double(get(hObject,'String')) returns contents of email as a double


% --- Executes during object creation, after setting all properties.
function email_CreateFcn(hObject, eventdata, handles)
% hObject    handle to email (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function piv_Callback(hObject, eventdata, handles)
% hObject    handle to piv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of piv as text
%        str2double(get(hObject,'String')) returns contents of piv as a double


% --- Executes during object creation, after setting all properties.
function piv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to piv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function institution_Callback(hObject, eventdata, handles)
% hObject    handle to institution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of institution as text
%        str2double(get(hObject,'String')) returns contents of institution as a double


% --- Executes during object creation, after setting all properties.
function institution_CreateFcn(hObject, eventdata, handles)
% hObject    handle to institution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function building_Callback(hObject, eventdata, handles)
% hObject    handle to building (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of building as text
%        str2double(get(hObject,'String')) returns contents of building as a double


% --- Executes during object creation, after setting all properties.
function building_CreateFcn(hObject, eventdata, handles)
% hObject    handle to building (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function room_Callback(hObject, eventdata, handles)
% hObject    handle to room (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of room as text
%        str2double(get(hObject,'String')) returns contents of room as a double


% --- Executes during object creation, after setting all properties.
function room_CreateFcn(hObject, eventdata, handles)
% hObject    handle to room (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cyverseUserName_Callback(hObject, eventdata, handles)
% hObject    handle to cyverseUserName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cyverseUserName as text
%        str2double(get(hObject,'String')) returns contents of cyverseUserName as a double


% --- Executes during object creation, after setting all properties.
function cyverseUserName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cyverseUserName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function password_Callback(hObject, eventdata, handles)
% hObject    handle to password (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of password as text
%        str2double(get(hObject,'String')) returns contents of password as a double


% --- Executes during object creation, after setting all properties.
function password_CreateFcn(hObject, eventdata, handles)
% hObject    handle to password (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on password and none of its controls.
function password_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to password (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if strcmp(eventdata.Key,'backspace')
    if ~isempty(handles.cyversePW)
        handles.cyversePW(end) = [];
    end
elseif strcmp(eventdata.Key,'shift')
    % no action
    eventdata.Key
elseif strcmp(eventdata.Key,'enter')
    
else
    handles.cyversePW = [handles.cyversePW eventdata.Character];
end


hv = repmat('*',[1 numel(handles.cyversePW)]);
set(handles.password,'String',hv);
v = get(handles.password,'String');
while ~strcmp(v,hv)
    pause(eps)
    set(handles.password,'String',hv);
    v = get(handles.password,'String');
end

guidata(hObject,handles)


% --- Executes on button press in register.
function register_Callback(hObject, eventdata, handles)
% hObject    handle to register (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get the cyverse user name
uname = get(handles.cyverseUserName,'String');
% get the cyverse password
pw = handles.cyversePW;
% init the cyverse account
cyverseConnect = myiinit_ver2(uname,pw);

% generate a random number
randN = num2str(randi(10^5),'%05d');
registrationData{1} = get(handles.name,'String');
registrationData{2} = get(handles.email,'String');
registrationData{3} = get(handles.piv,'String');
registrationData{4} = get(handles.institution,'String');
registrationData{5} = get(handles.building,'String');
registrationData{6} = get(handles.room,'String');
registrationData{7} = get(handles.cyverseUserName,'String');
registrationData{8} = randN;

% construct struct
rD = cell2struct(registrationData',{'Name','Email','PI','Institution','Building',...
    'Room','CyverseUname','randomN'});
uniqueID = [strrep(rD.Institution,' ','') '_' strrep(rD.PI,' ','') ...
    '_' strrep(rD.Building,' ','') '_' strrep(rD.Room,' ','') ...
    '_' strrep(rD.randomN,' ','')];
% add data to struct
rD.uid = uniqueID;
% add the path to the user name
rD.home = getenv('HOME');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if deployed then the below lines need to run
if isdeployed
    %CMD=['sed -i "s/uwtest0001/#uniqueID#/" /etc/condor/config.d/20-Name'];
    CMD=['sed -i "s/changeme/#uniqueID#/" /etc/condor/config.d/20-Name'];
    CMD = strrep(CMD,'#uniqueID#',rD.uid);
    system(CMD);


    %CMD=['sed -i "s/TEST/#randomN#/" /etc/condor/config.d/40-CCB'];
    %CMD = strrep(CMD,'#randomN#',rD.randomN);
    %system(CMD);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% get the file with ticket information
[file] = getTmpWriteTicket();
[~,nm] = fileparts(file);
tmp = ['/tmp/'];
cmd = ['iget '  file ' ' tmp];
[~,r] = system(cmd);
ticket = fileread([tmp nm]);
ticket(end) = [];
% write the JSON to local disk
fout = tempname;
regJSON = jsonencode(rD);
fileID = fopen(fout,'w');
fprintf(fileID,'%s',regJSON);
fclose(fileID);
[~,remoteNm] = fileparts(fout);
% push the registration data to remote service
remoteReg = ['/iplant/home/nmiller/prs/newFormRequest/'];
cmd = ['iput -t ' ticket ' ' fout ' ' remoteReg remoteNm];
[~,r] = system(cmd);
% wait for response
remoteReq = ['/iplant/home/nmiller/prs/newFormResponse/'];
cmd = ['echo -n ''' regJSON ''' | sha256sum '];
[~,uploadHash] = system(cmd);
uploadHash=uploadHash(1:end-4);
localResponseFile = [getenv('PHYTO_TMP') filesep uploadHash];
cmd = ['iget -f ' remoteReq uploadHash ' ' getenv('PHYTO_TMP') filesep];
[~,res] = system(cmd);
key = fileread(localResponseFile);
cmd = ['importKeys ' key];
[~,res] = system(cmd);

% --- Executes on key press with focus on register and none of its controls.
function register_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to register (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function register_CreateFcn(hObject, eventdata, handles)
% hObject    handle to register (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
