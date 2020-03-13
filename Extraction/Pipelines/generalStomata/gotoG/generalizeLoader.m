function [dataout] = generalizeLoader(fileName,loaderType,loaderArgs)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this is the generalized loader for a file the swtich statement determines the loader to be called
    % note that i am added functionality to this as follows:
    % if no loaderType is specified then i will look at the extension
    % this loader should be used with other piplines to load data into a
    % pipline and should be able to replace all loading funcnionality
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % infer loader type
    if nargin == 1
        fprintf(['Generalized loader is inferring datatype from extension\n']);
        imageTypeExtensions = {'nms','jpg','jpeg','tif','tiff','nef'};
        matTypeExtensionrs = {'mat'};
        [~,~,ext] = fileparts(fileName);
        ext(1) = [];
        if any(strcmp(ext,imageTypeExtensions))
            loaderType = 'simpleImageLoader';
        end
        if any(strcmp(ext,matTypeExtensionrs))
            loaderType = 'simpleMatLoader';
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%
    switch loaderType
        case 'stomata_histonormalize'
            dataout = imread(fileName);
            dataout = imhistmatch(dataout,loaderArgs{1});
        case 'simpleImageLoader'
            dataout = imread(fileName);
        case 'histogram'
            dataout = imread(fileName);
            dataout = imhistmatch(dataout,loaderArgs{1});
        case 'simpleMatLoader'
           % if isIRODS(fileName);
            dataout = load(fileName);
        case 'RGBsimpleImageLoader'
            dataout = imread(fileName);
            % remove the "alpha-channel"
            if size(dataout,3) > 3;dataout(:,:,4:end) = [];end
            % convert to double
            if isa(dataout,'uint8');dataout = double(dataout)/255;end
        case 'bvRGBsimpleImageLoader'
            dataout = imread(fileName);
            % remove the "alpha-channel"
            if size(dataout,3) > 3;dataout(:,:,4:end) = [];end
            % convert to double
            if isa(dataout,'uint8');dataout = double(dataout)/255;end
            % convert to double
            if isa(dataout,'uint16');dataout = double(dataout)/(2^16-1);end
    end
    % convert to double
    if isa(dataout,'uint8');dataout = double(dataout);end
    % convert uint16 to double
    if isa(dataout,'uint16');dataout = double(dataout);end
end
