function [] = plantProfiler(imageFile,debug,configFile)
    fprintf(['**************************************************\n']);
    fprintf(['*                Plant Profiler                  *\n']);
    fprintf(['*                  Version 1.02                 *\n']);
    fprintf(['*           Date Published: Jan 29, 2020         *\n']);
    fprintf(['*              Running at Biotron:works          *\n']);
    fprintf(['**************************************************\n']);
    if nargin <= 1;debug = false;end
    if nargin == 0;imageFile = [];end
    if isempty(imageFile);productionMode = true;else;productionMode = false;end
   
    
    qrCropBox = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['**************************************************\n']);
    CMD = 'echo $HOME';
    [status,home] = system(CMD);
    home(end) = [];
    fprintf(['Home location found as:' home '.\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if productionMode
        cameraSetup();
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % use the default config file location
    fprintf(['**************************************************\n']);
    if nargin <= 2
        % make default values
        configData{1} = Inf;
        configData{2} = Inf;
        configData{3} = Inf;
        configData{4} = Inf;
        configData{5} = Inf;
        configData{6} = Inf;
        configData{7} = Inf;
        configData{8} = 1;
        profilerScreenHome = '~/phytoMorphTK/plantProfiler/';
        mmkdir(profilerScreenHome);
        configFile = [profilerScreenHome 'default_config.csv'];
        
        
        if ~exist(configFile,'File')
            % write out default configue file for user
            cell2csv(configFile,configData);
        end
        
        
        fprintf(['Using default config file @:' configFile '\n']);
   
       
    end
    fprintf(['start: Reading config file @:' configFile '\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    configData = readtext(configFile);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read the config irods file for the CyVerse user
    fprintf(['end: Reading config file @:' configFile '\n']);
    fprintf(['**************************************************\n']);
    fprintf(['**************\n']);
    fprintf(['**************************************************\n']);
    fprintf(['start: Setting config parameters \n']);
    fprintf(['**************************************************\n']);
    %nm = '';
    fprintf(['Setting angle threshold of:' num2str(configData{1}) '\n']);
    angleThreshold = configData{1};
    fprintf(['Setting focal length target of:' num2str(configData{2}) '\n']);
    focalLengthTarget = configData{2};
    fprintf(['Setting focal length threshold of:' num2str(configData{3}) '\n']);
    focalLengthThreshold = configData{3};
    fprintf(['Setting f-number target of:' num2str(configData{4}) '\n']);
    fNumberTarget = configData{4};
    fprintf(['Setting f-number threshold of:' num2str(configData{5}) '\n']);
    fNumberThreshold = configData{5};
    fprintf(['Setting exposure-time target of:' num2str(configData{6}) '\n']);
    exposureTimeTarget = configData{6};
    fprintf(['Setting exposure-time  threshold of:' num2str(configData{7}) '\n']);
    exposureTimeThreshold = configData{7};
    fprintf(['Setting QR-search resize parameter of:' num2str(configData{8}) '\n']);
    re_resizeParameter = configData{8};
    fprintf(['end: Setting config parameters \n']);
    fprintf(['**************************************************\n']);
    fprintf(['**************\n']);
    fprintf(['**************************************************\n']);
    fprintf(['start: Reading CyVerse environment data for streaming service.\n']);
    fprintf(['**************************************************\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read the config irods file for the CyVerse user
    irodsJSONfile = '~/.irods/irods_environment.json';
    if ~exist(irodsJSONfile,'file')
        fprintf(['No cyverse user found. Please login.\n']);
        cmdLineLogin();
    end
    irodsJSON = loadjson('~/.irods/irods_environment.json');
    iPlantUser = irodsJSON.irods_user_name;
    fprintf(['found user: ' iPlantUser '.\n']);
    use = input(['Do you want to use this user name? (y,n)'],'s');
    if contains(use,'n');cmdLineLogin();end
    fprintf(['end: Reading CyVerse environment data for streaming service.\n']);
    fprintf(['**************************************************\n']);
    fprintf(['**************\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    try
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        basePicturePath = [home '/imageStore/'];
        fprintf(['**************************************************\n']);
        fprintf(['Setting up local image store @ ' basePicturePath '\n']);
        tmpLocalPath = [basePicturePath date filesep];
        mmkdir(tmpLocalPath);
        fprintf(['**************************************************\n']);     
        fprintf(['**************\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        remotePicturePath = '/iplant/home/#user#/maizeData/seedlingData/';
        remotePicturePath = strrep(remotePicturePath,'#user#',iPlantUser);
        tmpRemotePath = [remotePicturePath date filesep];
        CMD = ['imkdir -p ' tmpRemotePath];
        system(CMD);
        fprintf(['**************************************************\n']);
        fprintf(['Setting up remote image store base @: ' remotePicturePath '\n']);
        fprintf(['Setting up remote image store for session @: ' '\n']);
        fprintf(['**************************************************\n']);     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
       
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isdeployed
            fprintf(['**************************************************\n']);  
            fprintf(['start: Running in deployment mode and setting up JAR files \n']);
            fprintf(['**************************************************\n']);  
            file1 = [home '/phytoMorphTK/JARs/core-3.2.1.jar'];
            javaaddpath(file1);
            file2 = [home '/phytoMorphTK/JARs/javase-3.2.1.jar'];
            javaaddpath(file2);
            fprintf(['**************************************************\n']);  
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isempty(imageFile)
            fprintf(['**************************************************\n']); 
            imageFile = [home '/tmpData/tmpImage.nef'];
            [iPath,iName] = fileparts(imageFile);
            mmkdir(iPath)
            iPath = [iPath filesep];
            fprintf(['Setting temp image file name to: ' imageFile '\n']);
            fprintf(['**************************************************\n']); 
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        phytoStream = productionMode;
        if phytoStream
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % create datastream for: 
            %   "maizeSeedling" (1-input)" ; 
            %   "normal" (2-input)" 
            %   "1" - data stream 1
            %   "1" - active stream
            %   "@" the remote path
            cmd = 'clearAllStreams';
            system(cmd);
            fprintf(['**************************************************\n']);
            fprintf(['Setting up image stream.\n']);
            fprintf(['**************************************************\n']);
            createStream('plantProfiler','cyverseStream',1,tmpRemotePath,tmpLocalPath);
            %CMD = ['createStream maizeSeedling normal 1 1 "' tmpRemotePath '"'];
            %system(CMD);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % prompt user
        str = input('Press Enter for Picture and q for quit','s');
       
        
        Original_imageFile = imageFile;
        while isempty(str)
            totalTM = clock;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % call capture from phytoPhoto
            if productionMode
                captureImage('camera',1,iPath,iName);
                imageFile = fixNEFcase(Original_imageFile);
            end


            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf(['starting: image load, gray, and edge \n']);
            fprintf(['working on: ' imageFile '\n']);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % start processing the image
            if productionMode
                newName = [filesep home filesep '/tmpData/tmpImageTMP.nef'];
                CMD = ['cp "' imageFile '" "' newName '"'];
                system(CMD);
            else
                newName = imageFile;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % extract the meta data about the photo
            [p,nm] = fileparts(imageFile);
            info = imfinfo(imageFile);
            focalLength = info.DigitalCamera.FocalLength;
            fNumber = info.DigitalCamera.FNumber;
            exposureTime = info.DigitalCamera.ExposureTime;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % process for 
            % read the image
            I = imread(newName);
            % rectifiy
            [I,angle,nm,qrCropBox] = rectifyImage_ver3(I,re_resizeParameter,qrCropBox,debug);
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % print qr code to screen
            fprintf(['*******************************\n']);
            checkPass(1) = abs(angle) < angleThreshold;
            fprintf(['Angle:' num2str(angle) ':Check_Result:' num2str(checkPass(1)) ':Target:' num2str(angleThreshold) '\n']);
            fprintf(['*******************************\n']);
            checkPass(2) = abs(focalLength - focalLengthTarget) <= focalLengthThreshold;
            fprintf(['focalLength:' num2str(focalLength)  ':Check_Result:' num2str(checkPass(2)) ':Target:' num2str(focalLengthTarget) '\n']);
            fprintf(['*******************************\n']);
            checkPass(3) = abs(fNumber - fNumberTarget) <= fNumberThreshold;
            fprintf(['fNumber:' num2str(fNumber) ':Check_Result:' num2str(checkPass(3)) ':Target:' num2str(fNumberTarget) '\n']);
            fprintf(['*******************************\n']);
            checkPass(4) = abs(exposureTime - exposureTimeTarget) <= exposureTimeThreshold;
            fprintf(['exposureTime:' num2str(exposureTime) ':Check_Result:' num2str(checkPass(4)) ':Target:' num2str(exposureTimeTarget) '\n']);
            fprintf(['*******************************\n']);
            checkPass(5) = ~isempty(nm);
            fprintf(['QR Text:' nm '\n']);
            fprintf(['*******************************\n']);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if all(checkPass)
                if ~isempty(nm)
                    if productionMode
                        
                        % create the local base-date path if needed
                        newPath = [basePicturePath date filesep];
                        mmkdir(newPath)
                        newName = [newPath nm '.nef'];
                        CMD = ['cp "' imageFile '" "' newName '"'];
                        fprintf(['running command:' CMD '\n']);
                        system(CMD);
                        
                        
                        fprintf(['*******************************\n']);
                        fprintf(['Streaming file.\n']);
                        CMD = ['streamFile plantProfiler cyverseStream 1 "' newName '"'];
                        system(CMD);
                        fprintf(['*******************************\n']);
                        
                       
                        
                        
                        %{
                        % this is replaced by streams on Feb 18, 2019
                        tmpRemotePath = [remotePicturePath date filesep];
                        fprintf(['running command:' CMD '\n']);
                        CMD = ['imkdir -p "' tmpRemotePath '"'];
                        system(CMD);
                        %}
                        %{
                        % this is replaced by streams on Feb 19,2019
                        remoteName =  [tmpRemotePath nm '.nef'];
                        fprintf(['running command:' CMD '\n']);
                        CMD = ['iput -f "' newName '" "' remoteName '"'];
                        system(CMD);
                        %}

                        
                      
                    end
                    
                end
            else
                fprintf(['FAILED checks. Please retake picture.\n']);
            end
            
            totalTM = etime(clock,totalTM);
            fprintf(['Total Time: ' num2str(totalTM) '\n'])
            str = input('Press Enter for Picture and q for quit','s');
        end
    catch ME
        getReport(ME)
    end
    
    
end

%{
    testImage = '/mnt/scratch1/JUNK/plot256.nef';
    %testImage = '/home/nate/Downloads/tmpImage.nef';
    configFile = '/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/Maize/maizeSeedling/GRAIN_config.csv';
    
    testImage = '/home/nate/JUNK/{Plot_894}{Experiment_44}{Planted_4-24-2017}{SeedSource_DI 2114-3}{SeedYear_2016}{Genotype_Mo17}{Treatment_Control}{PictureDay_16}.nef';
    testImage = '/home/nate/tmpImage.NEF';
    plantProfiler(testImage);

    mcc -m blueScreenCapture.m -d /mnt/scratch1/phytomorph_dev/Aquire/blueScreen
%}