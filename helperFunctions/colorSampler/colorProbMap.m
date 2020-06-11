function [I,Io] = colorProbMap(I,density)
    MX = 255;
    if ischar(I);I = imread(I);end
    
    % store the original image for return
    Io = I;
    % size of image
    szI = size(I);
    
    
    % reshape samples by channels
    I = reshape(I,[prod(szI(1:2)) szI(3)]);
    I = density.binMapper(I);
    
    I = reshape(I,szI(1:2));
end