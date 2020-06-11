function [scores,error,errorL,U,E,L] = sampleAndDecompose(FileList,imageScale,funcSample,decomposeDIM)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read and sample the images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['Start:Loading images for creating the PCA.\n']);
    errorRM = [];
    % load first frame and pre-allocate data block
    tmp = imread(FileList{1});
    oSZ = size(tmp);
    tmp = imresize(tmp,imageScale);
    [tmpColorSample] = funcSample(tmp);
    I = zeros(numel(tmpColorSample),numel(FileList));
    totalToSample = numel(FileList);
    parfor e = 1:totalToSample
        try
            tmp = imread(FileList{e});
            oSZ = size(tmp);
            tmp = imresize(tmp,imageScale);
            [tmpColorSample] = funcSample(tmp);
            I(:,e) = tmpColorSample;
            fprintf(['Done reading:' num2str(e) ':' num2str(totalToSample) '\n']);
            rm(e) = false;
        catch ME
            rm(e) = true;
            getReport(ME)
        end
    end
    errorRM = find(rm);
    % remove the error images
    I(:,errorRM) = [];
    fprintf(['End:Loading images for creating the PCA.\n'])
    % decompose the samped image data
    [U,E,L] = PCA_FIT_FULL_Tws(I,decomposeDIM);
    scores = PCA_REPROJ_T(I,E,U);
    sims = PCA_BKPROJ_T(scores,E,U);
    error = sum((I - sims).*(I - sims),1);
    errorL = mean(error);
    %errorS = std(sum((I - sims).*(I - sims),1),1,2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end