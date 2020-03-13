function [data] = inLineFunctionExtractor_30102019(toOp,pointSet,opDomain,domainSize,extractionType,extractionArgs,prepareDislayImageParameters,disp)
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % 10.30.2019. ramI has been altred.  rawI is now calculated from toOp image
    % a prepare image function is called with the parameters in the second
    % to last variable
    %%%%%%%%%%%%%%%%%%%%
    % rawI := the raw image - for display purposes this image may need to
    % be altered to perpare it for overlay - for example the edges may need
    % to be trimmed off
    % trimValues := values to trim around the edges of the image - removed
    % and packaged  in second to last variable
    %%%%%%%%%%%%%%%%%%%%
    % toOp := image to operate on
    % pointSet := points to sample
    % 
   
    
    
    
    [rawI] = prepareImageforDisplay(toOp,prepareDislayImageParameters{1},prepareDislayImageParameters{2});
    
    disp = false;
    if disp
        figH = figure;
    else
        figH = '';
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the inLine function for loading data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set download options
    options = weboptions('Timeout',60);
    % set local web path for saving loader program
    webPath = './cachePath/';
    mkdir(webPath);
    % set loader program
    inlineFeatureExtractionProgram = [webPath 'generalizeFeatureExtractor.mat'];
    % save generalized loader from irods
    websave(inlineFeatureExtractionProgram,'https://de.cyverse.org/dl/d/5C41DB00-9619-4FB0-A18B-66757FC1503A/generalizeFeatureExtractor.mat',options);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load the loader
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % loader generalized loader
    featureP = load(inlineFeatureExtractionProgram);
    extractor = featureP.obj.func;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % peform extraction and analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['**********************************\n']);
    fprintf(['Starting opFunc over pointSet:\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for each slice
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %for slice = 1:size(toOp,3)
        
        %tmpOp = toOp(:,:,slice);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % for each point
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['start:Whole extracting bug eye image(s)\n']);
        str = clock;
        
        %%%%%%%%%%%%%%%%%%
        subI = generalizedSampler(toOp,pointSet(1,:),opDomain,domainSize,false);
        
        
        subI = thawTensor(subI,2);
        [T,NM] = extractor(subI,extractionType,extractionArgs);
        T = zeros(size(T,1),size(pointSet,1));
        
        %%%%%%%%%%%%%%%%%%
        parfor p = 1:size(pointSet,1)
            fprintf(['start:Extracting bug eye image ' num2str(p) ':' num2str(size(pointSet,1)) '.\n']);tic
            
            
            % start with sampler
            subI = bugEyeViaT(toOp,pointSet(p,:),opDomain,domainSize,disp,figH);
            subI = thawTensor(subI,2);
            
            
            % call extractor : a sequence of function layers
            [T(:,p)] = extractor(subI,extractionType,extractionArgs);
            
            
            eTime = toc;
            fprintf(['end:Extracting bug eye image ' num2str(p) ':' num2str(size(pointSet,1))  '@' num2str(eTime) '.\n']);
            
            
            rTime = mean(eTime)*(size(pointSet,1)-p);
            fprintf(['esti: time remaining : ' num2str(rTime) '\n']);
        end
        fprintf(['end:Whole extracting bug eye image(s) ' num2str(etime(clock,str)) '\n']);
        
        
            
        %T = opFunc(subI);
        %T = newF;
        %T = glueFrozenTensors(subI,newF);
        
    data.(NM) = T;
    data.rawI = rawI;
        
        
    fprintf('\n');
    fprintf(['Ending FFT:\n']);
    fprintf(['**********************************\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
