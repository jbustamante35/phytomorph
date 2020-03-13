function [I] = generateQRtile(str,sz)
    tmpFile = [tempname '.png'];
    CMD = ['qrencode ' str ' -o ' tmpFile];
    [r,o] = system(CMD);
    I = imread(tmpFile);
    %if size(I,1) < sz(1);I = imresize(I,sz,'method','nearest');end
    I = imresize(I,sz,'method','nearest');
    I = double(I);
    I = cat(3,I,I,I);
    delete(tmpFile);
end