function [] = captureMaizeSeedlingImages(imageFile,configFile)
    


    if nargin < 2
        configFile = '/phytomorph/GRAIN_config.csv';
    end

    configData = readtext(configFile);
    

    nm = '';
    angleThreshold = configData(1);
    focalLengthTarget = configData(2);
    focalLengthThreshold = configData(3);
    fNumberTarget = configData(4);
    fNumberThreshold = configData(5);
    exposureTimeTarget = configData(6);
    exposureTimeThreshold = configData(7);
    
    
    % read the config irods file for the CyVerse user
    irodsJSON = loadjson('~/.irods/irods_environment.json');
    iPlantUser = irodsJSON.irods_user_name;
    
    
    
    
    try
        
        
        basePicturePath = '/imageStore/';
        remotePicturePath = '/iplant/home/#user#/maizeData/seedlingData/';
        remotePicturePath = strrep(remotePicturePath,'#user#',iPlantUser);
        tmpRemotePath = [remotePicturePath date filesep];
        CMD = ['createStream maizeSeedling normal 1 1 "' tmpRemotePath '"'];
        system(CMD);
        
        
        fprintf(['RemotePath:' remotePicturePath '\n']);

        if isdeployed
            file1 = [filesep 'phytomorph' filesep 'core-3.2.1.jar']
            javaaddpath(file1);
            file2 = [filesep 'phytomorph' filesep  'javase-3.2.1.jar']
            javaaddpath(file2);
        end

        if isempty(imageFile)
            imageFile = [filesep home filesep '/tmpData/tmpImage.nef'];
            fprintf(['Testing on:' imageFile '\n']);
        end


        str = input('Press Enter for Picture and q for quit','s');

        while isempty(str)
            CMD = [filesep 'phytomorph' filesep 'MNsingleCaptureImage.sh'];
            fprintf(['Running Command:' CMD '\n\n']);
            [status,echoD] = system(CMD,'-echo');
            fprintf(['starting: image load, gray, and edge \n']);
            fprintf(['working on: ' imageFile '\n']);
            
            
            newName = [filesep home filesep '/tmpData/tmpImageTMP.nef'];
            CMD = ['cp "' imageFile '" "' newName '"'];
            system(CMD);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % path and name
            [p nm] = fileparts(imageFile);
            info = imfinfo(imageFile);
            focalLength = info.DigitalCamera.FocalLength;
            fNumber = info.DigitalCamera.FNumber;
            exposureTime = info.DigitalCamera.ExposureTime;
            % read the image
            I = imread(newName);
            % rectifiy
            [I,angle] = rectifyImage(I);
            try
                % get QR code
                nm = getQRcode_2(I);
            catch
                nm = [];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % print qr code to screen
            fprintf(['*******************************\n']);
            checkPass(1) = abs(angle) < angleThreshold;
            fprintf(['Angle:' num2str(angle) ':Check_Result:' num2str(checkPass(1)) ':Target' num2str(angleThreshold) '\n']);
            fprintf(['*******************************\n']);
            checkPass(2) = abs(focalLength - focalLengthTarget) <= focalLengthThreshold;
            fprintf(['focalLength:' num2str(focalLength)  ':Check_Result:' num2str(checkPass(2)) ':Target' num2str(focalLengthTarget) '\n']);
            fprintf(['*******************************\n']);
            checkPass(3) = abs(fNumber - fNumberTarget) <= fNumberThreshold;
            fprintf(['fNumber:' num2str(fNumber) ':Check_Result:' num2str(checkPass(3)) ':Target' num2str(fNumberTarget) '\n']);
            fprintf(['*******************************\n']);
            checkPass(4) = abs(exposureTime - exposureTimeTarget) <= exposureTimeThreshold;
            fprintf(['exposureTime:' num2str(exposureTime) ':Check_Result:' num2str(checkPass(4)) ':Target' num2str(exposureTimeTarget) '\n']);
            fprintf(['*******************************\n']);
            checkPass(5) = ~isempty(nm);
            fprintf(['QR Text:' nm '\n']);
            fprintf(['*******************************\n']);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if all(checkPass)
                if ~isempty(nm)
                    newPath = [filesep home basePicturePath date filesep];
                    CMD = ['mkdir -p "' newPath '"'];
                    fprintf(['running command:' CMD '\n']);
                    system(CMD);
                    newName = [newPath nm '.nef'];
                    CMD = ['cp "' imageFile '" "' newName '"'];
                    fprintf(['running command:' CMD '\n']);
                    system(CMD);
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
                    CMD = ['streamFile maizeSeedling normal 1 "' newName '"'];
                    system(CMD);
                    
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
    testImage = '/mnt/scratch1/JUNK/plot256.nef';
    configFile = '/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/Maize/maizeSeedling/GRAIN_config.csv';
    quickCheck_ver2(testImage,configFile);

    
    mcc -m captureMaizeSeedlingImages.m -d /mnt/scratch1/phytomorph_dev/Deploy/quickCheckMaizeSeedlings/
%}