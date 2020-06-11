function [sample] = sampleNhoodFeatures(I,samX)


    % read image
    I = imread(I);
    
    % get mean
    uI = imfilter(I,fspecial('disk',21));
    
    
    % size image
    %szI = size(uI);
    % reshape data
    %uI = reshape(uI,[prod(szI(1:2)) szI(3)]);
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % where to sample
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1: generate
    if isnumeric(samX)
        % generate random sample
        location = randperm(size(I,1));
        % sample 
        if N < 1;N = round(N*numel(location));end
        % pin down location
        location = location(1:N);
    % 2: load
    else
        location = samX.location;
    end
    
    % sample using the other function
    [sample] = sampleColorFeatures1(uI,samX);
    
end