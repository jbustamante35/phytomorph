function [] = set_and_captureImage(port,imageName,cameraSettings)
    setCameraParameters('f-number',cameraSettings(1));
    setCameraParameters('shutterspeed',cameraSettings(2));
    [capturePath,captureName] = fileparts(imageName);
    CMD = ['captureSingleAtPort_camera ' num2str(port) ' '''  capturePath ''' ''' captureName ''];
    CMD
end

%{
    set_and_captureImage(1,'~/testImage.tif',[3 4]);
%}