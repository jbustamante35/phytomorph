function [dataout] = generalizeLoader_301019(fileName,loaderType,loaderArgs)
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
        [~,~,ext] = fileparts(fileName);
        if any(strcmp(ext,imageTypeExtensions))
            loaderType = 'simpleImageLoader';
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%
    switch loaderType
        case 'stomata_histonormalize'
            dataout = imread(fileName);
            dataout = imhistmatch(dataout,loaderArgs{1});
        case 'simpleImageLoader'
            dataout = imread(fileName);
    end
end
%{
    % OLD LOADER FOR BACKUP
    function [I] = generalizeLoader(fileName,loaderType,loaderArgs)
    switch loaderType
        case 'stomata_histonormalize'
            I = imread(fileName);
            I = imhistmatch(I,loaderArgs{1});
    end
    end
%}
