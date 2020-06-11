function [fI] = testFilter(fileName)
    I = imread(fileName);
    fI = imfilter(I,fspecial('average',[10 10]));
end