function [STACK] = czi_tif_loader(fileName,filt)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PURPOSE:
    % load all frames from either a CZI 
    % or TIF stack 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUTS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fileName      := name of the stack file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the extension
    [fn,pth,ext] = fileparts(fileName);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if the file is a CZI file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(ext,'.czi')
        % open the file
        fprintf('Opening CZI stack.')
        M = bfopen(fileName);
        % STACK data
        cnt = 1;
        % make the STACK variable
        STACK = zeros([size(M{1}{1,1}) size(M{1},1)],'single');
        % loop over the frames
        for t = 1:1:size(M{1},1)
            % get the t-th frame
            tmp = single(M{1}{t,1});
            % THIS IS A FIX FOR THE DATA
            % set some of the max values to the min of the data frame
            tmp(tmp == max(tmp(:))) = min(tmp(:));
            % filter the data
            if ~isempty(filt);tmp = imfilter(tmp,filt,'replicate');end
            % log the data into STACK variable
            STACK(:,:,cnt) = tmp;
            % increase the cnt var
            cnt = cnt + 1;
            fprintf(['Done loading frame ' num2str(t) ':' num2str(size(M{1},1)) '\n']);
        end

        % convert to single
        STACK = single(STACK)/255;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if the stack is a TIF stack
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        fprintf('Opening TIF stack.')
        info = imfinfo(fileName);
        % for each frame
        for e = 1:numel(info)
            % get the MAX_N
            if e == 1
                tmp = imread(fileName,e);
                if isa(tmp,'uint16')
                    MAX_N = 2^16-1;
                else
                    MAX_N = 2^8-1;
                end
            end
            % read the data
            tmp = single(imread(fileName,e));
             % fix for the color data
            if size(tmp,3) > 1
                tmp = rgb2gray(tmp/255);
            end
            % THIS IS A FIX FOR THE DATA
            % set some of the max values to the min of the data frame
            tmp(tmp == max(tmp(:))) = min(tmp(:));
            % filter the data
            if ~isempty(filt);tmp = imfilter(tmp,filt,'replicate');end
            % increase the cnt var
            STACK(:,:,e) = tmp;
        end
        MAX_N = max(STACK(:));
        % convert to single
        STACK = single(STACK)/MAX_N;
    end
end