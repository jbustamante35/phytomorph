function [] = initNewDevice()
    % request user to unplu devices
    fprintf(['No devices are currently configured on device.\n']);
    input(['Please unplug all devices and hit enter.\n'],'s');
    % query the devices plugged in
    CMD = ['matlab_usbPortScanner1'];
    [r,o] = system(CMD);
    % ask the user to turn on device
    input(['Please plug in and turn on device. Hit enter when finished.\n'],'s');
    % scan usb ports
    CMD = ['matlab_usbPortScanner2'];
    [status,result] = system(CMD);
    % read in the new port data
    portData = readtext('~/tmpNewPort');
    portID = portData{1};
    % create the port via CLI 
    createPort(portID);
    fprintf(['Device added @ port:' portID '\n']);
end