function [sample] = sampleFeatures1(data,N)
    % default samN
    if nargin < 2;samN = .25;end
    % read image
    I = imread(data);
    % size image
    szI = size(I);
    % reshape data
    I = reshape(I,[prod(szI(1:2)) szI(3)]);
    % generate random sample
    ridx = randperm(size(I,1));
    % sample 
    if N < 1;N = round(N*numel(ridx));end
    % sample
    rgb = I(ridx(1:N),:);
    % generate hsv data
    hsv = rgb2hsv(rgb);
    % generate lab data
    lab = rgb2lab(rgb);
    % return features
    sample = [rgb,hsv,lab];
end