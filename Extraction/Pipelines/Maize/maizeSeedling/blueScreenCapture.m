function [] = blueScreenCapture(imageFile,debug,configFile)
    if nargin <= 1;debug = false;end
    fprintf(['**************************************************\n']);
    cameraSetup()
    fprintf(['**************************************************\n']);
    
    
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
        blueScreenHome = '~/phytoMorphTK/blueScreen/';
        mmkdir(blueScreenHome);
        configFile = [blueScreenHome 'default_config.csv'];
        
        
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
    irodsJSON = loadjson('~/.irods/irods_environment.json');
    iPlantUser = irodsJSON.irods_user_name;
    fprintf(['found user: ' iPlantUser '.\n']);
    fprintf(['end: Reading CyVerse environment data for streaming service.\n']);
    fprintf(['**************************************************\n']);
    fprintf(['**************\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    try
        
        
        
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        basePicturePath = '~/imageStore/';
        fprintf(['**************************************************\n']);
        fprintf(['Setting up local image store @ ' basePicturePath '\n']);
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
        fprintf(['Setting up remote image store for session @: ' tmpRemotePath '\n']);
        fprintf(['**************************************************\n']);     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        phytoStream = true;
        if phytoStream
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % create datastream for: 
            %   "maizeSeedling" (1-input)" ; 
            %   "normal" (2-input)" 
            %   "1" - data stream 1
            %   "1" - active stream
            %   "@" the remote path
            fprintf(['**************************************************\n']);
            fprintf(['Setting up image stream.\n']);
            fprintf(['**************************************************\n']);
            createStream('plantProfiler','cyverseStream',1,tmpRemotePath);
            %CMD = ['createStream maizeSeedling normal 1 1 "' tmpRemotePath '"'];
            %system(CMD);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['**************************************************\n']);
        CMD = 'echo $HOME';
        [status,home] = system(CMD);
        home(1) = [];
        home(end) = [];
        fprintf(['Home location found as:' home '.\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isdeployed
            fprintf(['**************************************************\n']);  
            fprintf(['start: Running in deployment mode and setting up JAR files \n']);
            fprintf(['**************************************************\n']);  
            file1 = [filesep 'phytomorph' filesep 'core-3.2.1.jar'];
            javaaddpath(file1);
            file2 = [filesep 'phytomorph' filesep  'javase-3.2.1.jar'];
            javaaddpath(file2);
            fprintf(['**************************************************\n']);  
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isempty(imageFile)
            fprintf(['**************************************************\n']); 
            productionMode = true;
            imageFile = [filesep home filesep '/tmpData/tmpImage.nef'];
            [iPath,iName] = fileparts(imageFile);
            iPath = [iPath filesep];
            fprintf(['Setting temp image file name to: ' imageFile '\n']);
            fprintf(['**************************************************\n']); 
        else
            productionMode = false;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % prompt user
        str = input('Press Enter for Picture and q for quit','s');

        while isempty(str)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % call capture from phytoPhoto
            if productionMode
                captureImage('camera',1,iPath,iName);
            end


            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf(['starting: image load, gray, and edge \n']);
            fprintf(['working on: ' imageFile '\n']);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % start processing the image
            if ~productionMode
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
            [I,angle,nm] = rectifyImage_ver2(I,re_resizeParameter,debug);
            
            
            
            
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
                    if ~productionMode
                        
                        
                        
                        %{
                        newPath = [filesep home basePicturePath date filesep];
                        CMD = ['mkdir -p "' newPath '"'];
                        fprintf(['running command:' CMD '\n']);
                        system(CMD);
                        newName = [newPath nm '.nef'];
                        CMD = ['cp "' imageFile '" "' newName '"'];
                        fprintf(['running command:' CMD '\n']);
                        system(CMD);
                        %}
                        
                        
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

                        
                        %CMD = ['streamFile maizeSeedling 1 normal "' newName '"'];
                        %system(CMD);
                        newName
                    end
                    
                end
            else
                fprintf(['FAILED checks. Please retake picture.\n']);
            end
            str = input('Press Enter for Picture and q for quit','s');
        end
    catch ME
        getReport(ME)
    end
    
    
end

%{

    


    mcc -m blueScreenCapture.m -d /mnt/scratch1/phytomorph_dev/Aquire/blueScreen
%}