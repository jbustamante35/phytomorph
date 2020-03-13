function [] = removeAllPorts()
    fprintf(['****************************************************\n']);
    fprintf(['Removing USB cameras from the system.\n']);
    fprintf(['****************************************************\n']);
    input('Make sure all your needed devices are plugged in and on! Then hit enter.','s');
    CMD = ['removePorts -silent'];
    system(CMD);
    fprintf(['Done removing all USB camera devices from system.\n']);
    fprintf(['****************************************************\n']);
end