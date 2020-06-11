function [] = makeImageMask(fileName,oPath,varargin)
    % get file parts
    [pth,nm,ext] = fileparts(fileName);
    % read image
    I = double(imread(fileName)) / 255;
    % make mask
    M = thresholdTasselImage(I, varargin{1});
    % get output image name
    outName = [oPath nm '.tif'];
    % save mask
    imwrite(M,outName);
end