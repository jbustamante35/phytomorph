function [H] = channelSampler(FileList,N_samFiles,N_samPoints,N_spaceVecs,imageLoader)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % input1: FileList
    % takes a list of files as input - this was a collection of image files
    % this is the default behavior
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % input: N_samFiles
    % the number of images/files to sample
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % input: N_samPoints


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % signature for image loader is:
    % (datafile,sampleNumber)
    % signature for mask is:
    % (datafile,sampleNumber,data);
    if nargin <= 4
        % default to using imread
        imageLoader{1} = @(x,ni)imread(x);
        % restrict the samples to a subset of the potential
        imageLoader{2} = @(x,ni,I)1:(size(I,1)*size(I,2));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % build out 8-bit space
    if nargin <= 3
        % assume 8-bit samples
        MX = 2^8-1;
        % read test file to gauage the size
        I = imageLoader{1}(FileList{1},1);
        % size of data
        szI = size(I);
        % number of channel
        numChannels = szI(end);
        % space vector samples
        for e = 1:numChannels
            N_spaceVecs{e} = linspace(0,MX,MX);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    % generate the size of the density space
    for e = 1:numel(N_spaceVecs)
        ds(e) = numel(N_spaceVecs{e});
    end
   
    % generate the output space
    density = zeros(ds);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % auto build the discretization and index mapper
    toBuild = true;
    if toBuild
        % apply discreteization to each sample
        df = @(X,i)discretize(X(:,i),N_spaceVecs{i});
        % loop over each dim
        g = @(x)arrayfun(@(i)df(x,i),1:numel(ds),'UniformOutput',false);
        % map the "rounded" sample the the linear-bin number
        binMapper = @(x)sub2ind(ds,cell2mat(g(x+1)));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
      
    % generate the random list
    ni = randi(numel(FileList),N_samFiles,1);
    
    
    for e = 1:numel(ni)
        try
            
            % load the image
            I = imageLoader{1}(FileList{ni(e)},ni(e));
            
            % load the potential sample points
            G = imageLoader{2}(FileList{ni(e)},ni(e),I);
            

            % size of image
            szI = size(I);
            % reshape into [samples x channels]
            I = reshape(I,[prod(szI(1:end-1)) szI(end)]);
            
            
            % generate random sample numbers
            % NOT: this is uniform over space
            sidx = randi(numel(G),N_samPoints,1);
            sidx = G(sidx);
            
            
            
            % record the percent of image sampled
            per = 100*numel(sidx)/numel(G);
            % sample
            I = I(sidx,:);
            % map sample in channel space to bin
            I = binMapper(I);
            % remove NaN
            I(isnan(I)) = [];
            % populate nd-histogram
            for e = 1:numel(I)
                density(I(e)) = density(I(e)) + 1;
            end
            % report
            fprintf(['sampled ' num2str(per) '\n']);
            
        catch ME
            getReport(ME)
        end
    end
    density = density / sum(density(:));
    H.density = density;
    H.axis = N_spaceVecs;
    H.binMapper = binMapper;
    H.prob = @(x)prob(x,H);
end