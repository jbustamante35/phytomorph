function [] = sampleStructureFeatures(I)
    % read image
    I = imread(I);
    % convet to gray
    if size(I,3) > 1;G = rgb2gray(I);end
    
    
end