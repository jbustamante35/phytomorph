function [sample] = sampleColorFeatures1(I,samX)
    % default samN
    if nargin < 2;samX = .25;end
    
    % read image
    I = imread(I);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % what to sample
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % size image
    szI = size(I);
    % reshape data
    I = reshape(I,[prod(szI(1:2)) szI(3)]);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % where to sample
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate
    if isnumeric(samX)
        % generate random sample
        location = randperm(size(I,1));
        % sample 
        if samX < 1;samX=round(samX*numel(location));end
        % pin down location
        location = location(1:samX)';
    % load
    else
        location = samX.location;
    end
    
    
    % sample
    rgb = I(location,:);
    % convert
    if isa(rgb,'uint8');rgb = double(rgb)/255;end
    % generate hsv data
    hsv = rgb2hsv(rgb);
    % generate lab data
    lab = rgb2lab(rgb);
    
    
    if isa(samX,'struct')
        % return features data
        samX.data = [samX.data, rgb,hsv,lab];
        sample = samX;
    else
        % return features data
        sample.data = [rgb,hsv,lab];
        % return location data
        sample.location = location;
    end
   
    
end