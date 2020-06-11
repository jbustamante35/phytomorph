function [I] = sampleImageList(FileList,imageScale,funcSample)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read and sample the images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['Start:Loading images.\n']);
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
    fprintf(['End:Loading images.\n']);
end