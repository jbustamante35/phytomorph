function [] = captureImage(type,usbPortId,oPath,oFile)
    fprintf(['Capturing image.\n']);
    CMD = ['captureSingleImage ' type ' ' num2str(usbPortId) ' ' oPath ' ' oFile];
    system(CMD);
end