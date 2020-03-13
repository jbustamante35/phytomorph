function varargout = deviceBank_GUI(varargin)
% DEVICEBANK_GUI MATLAB code for deviceBank_GUI.fig
%      DEVICEBANK_GUI, by itself, creates a new DEVICEBANK_GUI or raises the existing
%      singleton*.
%
%      H = DEVICEBANK_GUI returns the handle to a new DEVICEBANK_GUI or the handle to
%      the existing singleton*.
%
%      DEVICEBANK_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DEVICEBANK_GUI.M with the given input arguments.
%
%      DEVICEBANK_GUI('Property','Value',...) creates a new DEVICEBANK_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before deviceBank_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to deviceBank_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help deviceBank_GUI

% Last Modified by GUIDE v2.5 24-Jan-2020 19:01:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @deviceBank_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @deviceBank_GUI_OutputFcn, ...
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


% --- Executes just before deviceBank_GUI is made visible.
function deviceBank_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to deviceBank_GUI (see VARARGIN)

% Choose default command line output for deviceBank_GUI
handles.output = hObject;


home = getenv('HOME');
reg = csvread([home '/phytoMorphTK/reg.txt']);


if ~reg
    
    registrationData = inputdlg({'Name','Institution','Principal Investigator','Building','Room #','email'});
    [uname,pw] = logindlg('Title','Main cyverse acount');
    randN = num2str(randi(10^5),'%05d');
    registrationData{end+1} = uname;
    registrationData{end+1} = randN;
    myiinit_ver2(uname,pw);
    rD = cell2struct(registrationData,{'Name','Institution','PI','Building','Room','Email','CyverseUname','randomN'});
    uniqueID = [strrep(rD.Institution,' ','') '_' strrep(rD.PI,' ','') ...
        '_' strrep(rD.Building,' ','') '_' strrep(rD.Room,' ','') ...
        '_' strrep(rD.randomN,' ','')];
    rD.uid = uniqueID;
    
    
    CMD=['sed -i "s/uwtest0001/#uniqueID#/" /etc/condor/config.d/20-Name'];
    CMD = strrep(CMD,'#uniqueID#',rD.uid);
    system(CMD);
    
    
    CMD=['sed -i "s/TEST/#randomN#/" /etc/condor/config.d/40-CCB'];
    CMD = strrep(CMD,'#randomN#',rD.randomN);
    system(CMD);
    
  
    
    
    
end



% mode op
devMode = true;        

% init the string buffer
handles.buffer = stringBuffer();

% run the device setup block
if ~devMode
    handles.numberDevices = deviceSetup();
else
    handles.numberDevices = 2;
end

portList = 1:2;
handles.curSpark = struct;


% get screen resolution
handles.res = getScreenResolution();

% set the block size
handles.blockSZ = round(2*handles.res);

% init the current scan state
handles.curTile = 1;

% number of states per tile strip
handles.tilesPerStrip = 4;

% initState Matrix
iM = [[0;1;0],repmat([0;0;1],[1 handles.tilesPerStrip-1])];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% master color tile
% init the color states - hue will change along timing
% 1) end tiles 2) sample tiles 3) meta data tiles
initColorStates = [[66 173 84];[42 61 209];[167 66 245]]/255;
rgbStates = generateColorPanel(handles.tilesPerStrip,handles.numberDevices,initColorStates);
  
% init the saturation level - level for each state
initSat = [.5 1 .1];
satLevel = repmat(initSat,[size(rgbStates,1) 1]);



% init the state figure
machineState = ones(handles.numberDevices*handles.blockSZ(1),...
                handles.tilesPerStrip*handles.blockSZ(2));
machineState = cat(3,machineState,machineState,machineState);
for portIdx = 1:handles.numberDevices
    for stateIdx = 1:handles.tilesPerStrip
        satL = initSat*iM(:,stateIdx);
        machineState = updateStateImage(machineState,portIdx,stateIdx,...
        handles.blockSZ,rgbStates(stateIdx,:),satL,[],[]);
    end
end



% set the state to first scan
handles.firstScan = true;

% update the machine state image
handles.machineState = machineState;

% show the init state
axes(handles.axes1);
imshow(handles.machineState,[]);
   
% handles state image
handles.sImage = stateImage(machineState,handles.axes1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop over the devices to set up device state vectors
processQueue = [];
for e = 1:handles.numberDevices
    currentQueue = stateBlock(eye(2)',[],[],[],stateDirac(2,1));
    processQueue = [processQueue,currentQueue];
end

%%%%%%%%%%%%%%%%%%%%
%{
% sample data transistion condition
transitionFuncs{2} = @(S,X)containKVP(X,{'{dataType_sampleCode}'});
% meta data transistion condition
transitionFuncs{3} = @(S,X)containKVP(X,{'{dataType_metaCode}'});
% stop code transistion condition
transitionFuncs{4} = @(S,X)containKVP(X,{'{command_sessionEnd}'});
%}
% sample data transistion condition
transitionFuncsArgs{1} = {'{t_S}','{u_#}'};
% sample data transistion condition
transitionFuncsArgs{2} = {'{t_s}','{u_#}'};
% meta data transistion condition
transitionFuncsArgs{3} = {'{t_m}','{u_#}'};
% stop code transistion condition
transitionFuncsArgs{4} = {'{t_E}','{u_#}'};
%transFunc = @(curState,stateTarget,curMsg,msgTarget)all(curState==stateTarget)&containKVP(curMsg,msgTarget);
transFunc = @(curState,stateTarget,curMsg,msgTarget)b2q(containKVP(curMsg,msgTarget));
%%%%%%%%%%%%%%%%%%%%
targetState = [0;1;0];

successFunc{1} = @(S,X)sparkProcess(X,S,'begin');
successFunc{2} = @(S,X)sparkProcess(X,S,'sample');
successFunc{3} = @(S,X)sparkProcess(X,S,'metadata');
successFunc{4} = @(S,X)sparkProcess(X,S,'end');


failFuncs{1} = @(S,X)fprintf(['Wrong tile scan!']);
failFuncs{2} = @(S,X)fprintf(['Wrong tile scan!']);
failFuncs{3} = @(S,X)fprintf(['Wrong tile scan!']);
failFuncs{4} = @(S,X)fprintf(['Wrong tile scan!']);



%{
updateFunc{1} = updateStateImage(machineState,portIdx,stateIdx,...
        handles.blockSZ,handles.rgbStates(stateIdx,:),handles.satLevel(1),[],[]);
%}


handles.stateBlocks = [];
for portIdx = 1:handles.numberDevices
    
    %{
    tmpCondition = {'{command_sessionStart}','{usbPort_#}'};
    cmd = strrep(tmpCondition{2},'#',num2str(portList(portIdx)));
    % start command transition condition
    transitionFuncs{1} = @(S,X)containKVP(X,{cmd});
    %}
    
    tmpT = {};
    for e = 1:numel(transitionFuncsArgs)
        
        tmpArgs = transitionFuncsArgs{e};
        tmpArgs{2} = strrep(tmpArgs{2},'#',num2str(portList(portIdx)));
        % start command transition condition
        tmpT{e} = @(curState,curMsg)transFunc(curState,targetState,curMsg,tmpArgs);
    end
    
    
    
    curBlock = generateStateBlocks(handles.sImage,portIdx,...
        handles.tilesPerStrip,handles.numberDevices,...
        tmpT,successFunc,failFuncs);
    
    % attach a row-sequence of state blocks to string buffer
    handles.buffer.attachListener(curBlock);
    
    %for r = 1:numel(curBlock)
    %    handles.buffer.attachListener(curBlock(r));
    %end
    
    handles.stateBlocks = [handles.stateBlocks;curBlock];
end







% init states
handles.state = ones(handles.numberDevices,handles.tilesPerStrip);
        

buttonSZ = 201;
disk = zeros(buttonSZ);
hsz = hisize(disk);
disk(hsz(1),hsz(2)) = 1;
buttonRadius = 40;
msk = bwdist(disk) < buttonRadius;
CL = [[210,219,26];[97,176,7];[161,36,27]]/255;
buttons = .94*ones(buttonSZ,buttonSZ,3,size(CL,1));
cidx = find(msk);
for m = 1:size(CL,1)
    tmp = buttons(:,:,:,m);
    for k = 1:3
        tmp2 = tmp(:,:,k);
        tmp2(cidx) = CL(m,k);
        tmp(:,:,k) = tmp2;
    end
    buttons(:,:,:,m) = tmp;
end
handles.buttons = buttons;

axes(handles.axes2);
imshow(buttons(:,:,:,1),[]);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes deviceBank_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = deviceBank_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set focus to edit box
uicontrol(handles.edit1);

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on edit1 and none of its controls.
function edit1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
pause(eps);
if strcmpi(eventdata.Key,'return')
    
    
    
    % drain the buffer after filling
    str = eventdata.Source.String;
    while isempty(str)
        pause(eps);
        str = get(handles.edit1,'String');
    end
    % get the buffer before clearing it
    str = get(handles.edit1,'String');
    
    
    % read current buffer
    cur = get(handles.edit1,'String');
    while ~isempty(cur)
        % clear the buffer
        set(handles.edit1,'String','');
        % read buffer
        cur = get(handles.edit1,'String');
    end
   
    
    
    if validPSON(str)

        % increment the current tile on the strip
        handles.curTile = handles.curTile + 1;
        if handles.firstScan
            handles.firstScan = false;
        end

        axes(handles.axes2);
        imshow(handles.buttons(:,:,:,2),[]);
        
        pause(.5)
        axes(handles.axes2);
        imshow(handles.buttons(:,:,:,1),[]);
        
        % set focus to edit box
        uicontrol(handles.edit1);
        
        handles.buffer.string = str;
        
        % bubble process
        %handles.stateBlocks.processCommand(handles.curSpark,str)
    else
        
        axes(handles.axes2);
        imshow(handles.buttons(:,:,:,end),[]);
        
        pause(.5)
        axes(handles.axes2);
        imshow(handles.buttons(:,:,:,1),[]);
    end
   
end
guidata(hObject,handles);

%handles.buffer


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
answer = inputdlg('Please enter registration key','Register Computer');
