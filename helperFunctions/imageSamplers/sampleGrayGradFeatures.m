function [sample] = sampleGrayGradFeatures(I,sample)
    % read image
    I = imread(I);
    % convet to gray
    if size(I,3) > 1;G = rgb2gray(I);end
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % what to sample
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % calc grad
    if ~isa(G,'double');G = double(G);end
    [d1,d2] = gradient(G);
    grad = cat(3,d1,d2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % where to sample
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate
    if isnumeric(sample)
        % generate random sample
        location = randperm(size(I,1));
        % sample 
        if sample < 1;sample=round(sample*numel(location));end
        % pin down location
        location = location(1:sample)';
    % load
    else
        location = sample.location;
    end
    
    
    
    % sample using the other function
    [sample] = sampleColorFeatures1(grad,sample);
    
end