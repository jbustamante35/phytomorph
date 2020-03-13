function [] = setCameraParameters(parameter,value)
    CMD = ['gphoto2 --set-config ' parameter '=' num2str(value)];
    %system(CMD);
    CMD
end

%{
setCameraParameters('f-number','10')
    
%}