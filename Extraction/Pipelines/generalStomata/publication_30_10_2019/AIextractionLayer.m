function [data] = AIextractionLayer(fileName,loaderType,loaderArgs,extractionType,extractionArgs,domainData,pointSetData,customData,oPath,rPath)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extraction layer - file copied on 30.10.2019
    % This function will 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fileName := name of the file to operate on
    % loaderType := parameter for the loader
    % loaderArgs ;= loader parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    mkdir(oPath);
    
    pointSet = pointSetData;
    opDomain = domainData{1};
    domainSize = domainData{2};
    trimValues = domainData{3};
    
    
    prepareImageForDisplay = {'edgeTrim',trimValues};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % to mod the generalized loader
    % edit generalizeLoader.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % to publish a new version of the generalized loader
    % publishToCyVerse('gLoader');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the loader
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set download options
    options = weboptions('Timeout',60);
    % set local web path for saving loader program
    webPath = './cachePath/';
    mkdir(webPath);
    % set loader program
    loaderProgram = [webPath 'generalizeLoader.mat'];
    % save generalized loader from irods
    websave(loaderProgram,'https://de.cyverse.org/dl/d/B883590D-9C5F-470F-9738-C3948A6D44A4/generalizeLoader.mat',options);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load the loader
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % loader generalized loader
    loaderP = load(loaderProgram);
    loader = loaderP.obj.func;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % execute loader
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    I = loader(fileName,loaderType,loaderArgs);
    imageSZ = size(I);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    preP
    [data] = inLineFunctionExtractor_30102019(I,pointSet,opDomain,domainSize,extractionType,extractionArgs,prepareImageForDisplay,false);
    
    
    
    
    data.customData = customData;
    if ~isempty(oPath)
        fileList{1} = [oPath customData{1} '.mat'];
        save(fileList{1},'-struct','data');
        pushToiRods(rPath,fileList);
    end
    
end