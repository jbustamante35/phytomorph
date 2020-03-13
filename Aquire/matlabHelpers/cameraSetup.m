function [] = cameraSetup()
     
    str = input(['Do you want to reset all USB camera ports? (y/n)'],'s');
    if contains(lower(str),'y');removeAllPorts();end
    

    [o] = listPorts();
    if isempty(o)
        % set up the camera
        fprintf(['No cameras are currently configured on device.\n']);
        input(['Please unplug all devices and hit enter.\n'],'s');
        CMD = ['matlab_usbPortScanner1'];
        system(CMD);
        input(['Please plug in and turn on camera. Hit enter when finished.\n'],'s');
        CMD = ['matlab_usbPortScanner2'];
        [status,result] = system(CMD);
        portData = readtext('~/tmpNewPort');
        portID = portData{1};
        createPort(portID);
        fprintf(['Camera added @ port:' portID '\n']);
    else
        % checking for camera is not yet set up
        input(['Make sure camera is on and hit enter.'],'s');
    end
end


%{
    mcc -m cameraSetup.m -d /mnt/scratch1/phytomorph_dev/Aquire/blueScreen
%}