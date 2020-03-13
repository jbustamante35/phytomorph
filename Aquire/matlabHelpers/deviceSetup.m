function [n] = deviceSetup()
     
    % list the ports
    str = input(['Do you want to reset all USB device ports? (y/n)'],'s');
    if contains(lower(str),'y');removeAllPorts();end
    
    % list the ports
    [o] = listPorts();
    
     % if empty
    if isempty(o)
        initNewDevice();
    else
        % checking for camera is not yet set up
        input(['Make sure camera is on and hit enter.'],'s');
    end
    
    str = input(['Do you want to add another device? (y/n)'],'s');
    while ~contains(lower(str),'y')
        initNewDevice();
        str = input(['Do you want to add another device? (y/n)'],'s');
    end
    
    [o,n] = listPorts();
end


%{
    mcc -m cameraSetup.m -d /mnt/scratch1/phytomorph_dev/Aquire/blueScreen
%}