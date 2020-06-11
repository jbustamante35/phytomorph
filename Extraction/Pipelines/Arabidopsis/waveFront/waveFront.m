function [] = waveFront(fileName,oPath)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUTS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fileName      := name of the stack file
    % oPath         := location to write the results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % handle the file manipulations
    % make the output location
    mkdir(oPath)
    % get file parts
    [pth,nm,ext] = fileparts(fileName);
    % get mask images
    maskFolder = [pth filesep 'masks' filesep];
    FilePath = maskFolder;
    FileList = {};
    FileExt = {'tif'};
    % look for the mask files
    maskFileList = gdig(FilePath,FileList,FileExt,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set the configure data from file if present
    configData.globalSKIP = 1;
    % colors for displays
    configData.CL =  {'r' 'g' 'b'};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read the csv file
    if ~isdeployed()
        %configureFile = '/home/nate/WaveFrontConfigure.csv';
        configureFile = '';
    else
        configureFile = [pth filesep 'WaveFrontConfigure.csv'];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read the configure file from disk
    % unless it does not exist
    % if not exist then make one with default values
    if exist(configureFile)
        % read from the configure file
        fileData = readtext(configureFile);
        for col = 1:size(fileData,2)
            vecData = [];
            for row = 2:size(fileData,1)
                if ~isempty(fileData{row,col})
                    vecData = [vecData fileData{row,col}];
                end
            end
            configData.(fileData{1,col}) = vecData;
        end
    else
        % flag to render video to disk
        configData.toRenderVideo = 1; 
        % which are less than PER of the max
        configData.waveThresholdPercent = [.2 .35 .5];
        % skips for the frame over the t-stacks
        configData.t_stack_SKIP_frame_vec = 1:10:31;
        % percent contour levels to process the data
        configData.PER_contour = [.2 .35 .5];
        % filter default values
        configData.load_filt = [41 41 11];
        % filter for default time smoothing
        configData.time_filt = 11;
        % filter for default time smoothing
        configData.gradient_filt = [11 11 3];
        % load filter data
        configData.load_filt = [21 21 3];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load and stack data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % can be either a CZI of TIF file
    filt = fspecial('gaussian',configData.load_filt(1:2),configData.load_filt(3));
    STACK = czi_tif_loader(fileName,filt);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make whole plant mask
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make the std of the data across time
    SS = std(STACK,1,3);
    % make the mean signal of the image
    U = mean(STACK,3);
    % make the plant mask based on the standard dev
    plantMask = bindVec(SS) > graythresh(bindVec(SS))*.2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the locations of the plant mask
    midx = find(plantMask);
    % sample the mean image @ the plant mask
    uI = U(midx);
    % split the mean values into two groups
    meanLevel = graythresh(uI);
    % high mask
    hMask = uI > meanLevel;
    % low mask
    lMask = uI < meanLevel;
    % high group index
    highIndex = midx(find(uI > meanLevel));
    % low group index
    lowIndex = midx(find(uI < meanLevel));
    % make high mask
    highMask = zeros(size(plantMask));
    highMask(highIndex) = 1;
    % make low mask
    lowMask = zeros(size(plantMask));
    lowMask(lowIndex) = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write out the masks
    imwrite(double(255*SS),[oPath 'std_image.tif']);
    imwrite(double(255*U),[oPath 'mean_image.tif']);
    imwrite(double(plantMask),[oPath 'std_mask.tif']);
    imwrite(double(highMask),[oPath 'mean_highMask.tif']);
    imwrite(double(lowMask),[oPath 'mean_lowMask.tif']);
    imwrite(double(lowMask.*plantMask),[oPath 'mean_lowMask_AND_std_mask.tif']);
    imwrite(double(highMask.*plantMask),[oPath 'mean_highMask_AND_std_mask.tif']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % process each of the masks using the 'focus' from each of the maskFileList
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    innerLoopForMask(STACK,plantMask,'rawPlant',oPath,nm,maskFileList,configData);
    
    % changed to only calculate the raw mask
    %innerLoopForMask(STACK,highMask==1 & plantMask==1,'highMask',oPath,nm,maskFileList,configData);
    %innerLoopForMask(STACK,lowMask==1 & plantMask==1,'lowMask',oPath,nm,maskFileList,configData);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    here = 1;
    
    
    
end

%{

    FilePath = '/home/nate/Downloads/glr3.3/glr3.3/';
    FileList = {};
    FileExt = {'tif','TIF'};
    FileList = gdig(FilePath,FileList,FileExt,1);
    
    oPath = '/home/nate/Downloads/glr3.3/';
    oFile = [oPath 'glr3.3.tif'];

    for e = 1:numel(FileList)
        I = imread(FileList{e});
        imshow(I,[]);
        drawnow
        imwrite(I,oFile,'writemode','append');
    end
  
    waveFront(oFile,'./output/');





    FilePath = '/mnt/tetra/nate/arthur/waveFront/';
    FileList = {};
    FileExt = {'czi'};
    tic
    FileList = gdig(FilePath,FileList,FileExt,1);
    waveFront(FileList{1},'./output/');


    fileName = '/mnt/tetra/nate/arthur/Nobel Prize/649/Experiment-649.czi';
    fileName = '/mnt/snapper/nate/Test1/test1.tif';
    fileName = '/mnt/snapper/nate/ART_test/Exp149-Substack_49-621.tif';
    %fileName = '/mnt/snapper/nate/ART_test/Experiment-149.czi';
    waveFront(fileName,'./output/');
%}