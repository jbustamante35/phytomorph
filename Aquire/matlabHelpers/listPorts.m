function [o,n] = listPorts()
    CMD = ['listPorts'];
    [~,o] = system(CMD);
    fidx = strfind(o,char(10));
    n = numel(fidx);
end

%{



 mcc -m listPorts.m -d /mnt/scratch1/phytomorph_dev/Aquire/matlabHelpers

%}