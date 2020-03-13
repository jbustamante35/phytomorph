function [] = deviceBank(dev,debug)
    dev = true;
    debug = false;
    firstScan = true;
    close all
    try 
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['**************************************************\n']);
        fprintf(['*                Image Device Capture             *\n']);
        fprintf(['*                  Version 1.01                  *\n']);
        fprintf(['*           Date Published: Jan 13, 2020         *\n']);
        fprintf(['**************************************************\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['**************************************************\n']);
        CMD = 'echo $HOME';
        [status,home] = system(CMD);
        home(end) = [];
        fprintf(['Home location found as:' home '.\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        
        
        n=2;
        if ~dev;n = deviceSetup();end
        
        
        res = getScreenResolution();
        blockSZ = round(2*res);
        curTile = 1;
        statesPerTile = 4;
        initSampleNumber = 2;
        initColrStates = [[1 0 0];[0 1 0];[0 0 1]];
        [rgbStates] = generateColorPanel(statesPerTile,initSampleNumber,initColrStates);
       
        
        % init the saturation level
        initSat = [1 .5];
        satLevel = repmat(initSat,[size(rgbStates,1) 1]);
        
        
        % init states
        state = ones(n,statesPerTile);
        
       
        
        
        
        % init block size to 100 x 100
        statesPerTile = 4;
        machineState = ones(n*blockSZ(1),statesPerTile*blockSZ(2));
        machineState = cat(3,machineState,machineState,machineState);
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % NOTE: there is a problem if the last code is scanned before all the 
        % meta data tiles
        % NOTE: the sample codes can be scanned out of order if there are more
        % than 1
        %%%%%%%%%%%%%%%%%%%%%%%%
        initBase = ['imageStore' filesep];
        folderString =[getenv('HOME') filesep initBase];
        fileString = '';

        
        
        if nargin == 0;debug=false;end
        %%%%%%%%%%%%%%%%%%%%%%%%
        % init vars
        % init the prompts
        initPrompt = 'Please scan init command.\n';
        samplePrompt = 'Please scan sample code #.\n';
        metaPrompt = 'Please scan meta data tile number #.\n';
        endPrompt = 'Please scan end command.\n';

        %%%%%%%%%%%%%%%%%%%%%%%%
        % scan until the start code 
        % init the 
        str = '';
        %while ((~containsKeys(str,{'command','usbPort'})) || (~flag))
        while ~containKVP(str,{'{command_sessionStart}','{usbPort_*}'})
            str = getValidRead(initPrompt);
            if ~containKVP(str,{'{command_sessionStart}'})
                fprintf(['****************************\n']);
                fprintf(['Please scan the first tile!\n']);
                fprintf(['****************************\n']);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%
           
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % report the start code has been found
        sampleCodeN = str2num(getValue(findKVP(str,'sampleSpecific')));
        totalCodeN = str2num(getValue(findKVP(str,'totalTiles')));
        usbPort = getValue(findKVP(str,'usbPort'));
        deviceType = getValue(findKVP(str,'deviceType'));
        resolution = getValue(findKVP(str,'resolution'));
        fprintf(['Session started for ' deviceType '@' usbPort ' usbPort.\n']);
        expectedTilesN = totalCodeN;
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        
        totalVisN = totalCodeN + sampleCodeN + 2;
        tic;[qrImage] = generateQRtile(str,blockSZ);toc
        totalVisN = totalCodeN + sampleCodeN + 2;
        
        
        
        if totalVisN ~= statesPerTile
            statesPerTile = totalVisN;
            [rgbStates] = generateColorPanel(statesPerTile,initSampleNumber,initColrStates);
            machineState = ones(n*blockSZ(1),statesPerTile*blockSZ(2));
            machineState = cat(3,machineState,machineState,machineState);
        end
        
        %close all
        if firstScan
            for portIdx = 1:n
                for stateIdx = 1:totalVisN
                    machineState = updateStateImage(machineState,portIdx,stateIdx,blockSZ,rgbStates(stateIdx,:),satLevel(1));
                    %imshow(machineState,[]);
                    %drawnow
                end
            end
        end
       
        
        r = str2num(usbPort);
        machineState = updateStateImage(machineState,r,curTile,blockSZ,rgbStates(curTile,:),satLevel(1),qrImage);
        curTile = curTile + 1;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % scan sample codes
        %%%%%%%%%%%%%%%%%%%%%%%%
        sampleCnt = 1;
        sampleFileData = '';
        fprintf(['Expecting ' num2str(sampleCodeN) ' sample QR codes.\n'])
        while ((~containKVP(str,{'{dataType_sampleCode}'})) && (sampleCodeN ~= 0))
            tmpPrompt = strrep(samplePrompt,'#',num2str(sampleCnt));
            str = getValidRead(tmpPrompt);
            if containKVP(str,{'{dataType_sampleCode}'})
                % get file data from sample code
                sampleFileData = getValue(findKVP(str,'fe'));
                % if resolution should be in sample code
                if strcmp(resolution,'sc')
                    if containsKeys(str,'resolution')
                        resolution = getValue(findKVP(str,'resolution'));
                    else
                        fprintf(['Expected resolution in sample code.\n']);
                        fprintf(['Proceeding with default of 1200 dpi.\n']);
                        resolution = '1200';
                    end
                end
                % decrement
                sampleCodeN = sampleCodeN - 1;
                
                
                tic;[qrImage] = generateQRtile(str,blockSZ);toc
                machineState = updateStateImage(machineState,r,curTile,blockSZ,rgbStates(curTile,:),satLevel(1),qrImage);
                curTile = curTile + 1;
                
                
            else
                fprintf('Your command was not a valid sample QR.\n');
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%

        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % scan meta data codes until the top is scanned
        %%%%%%%%%%%%%%%%%%%%%%%%
        tileCount = 2;
        metaTileCount = 0;
        while (~(containKVP(str,{'{dataType_metaCode}'})) || (metaTileCount ~= expectedTilesN))
            tmpPrompt = strrep(metaPrompt,'#',num2str(tileCount));
            str = getValidRead(tmpPrompt);
            if containKVP(str,{'{dataType_metaCode}',['{tileNumber_' num2str(tileCount) '}']})
                % get the folder values
                folderData = getValue(findKVP(str,'fr'));
                % get the file values
                fileData = getValue(findKVP(str,'fe'));
                % if the proper tile is scanned
                folderString = [folderString folderData];
                fileString = [fileString fileData];
                % report if debug
                if debug
                    fprintf(['Folder:' folderString '\n']); 
                    fprintf(['File:' fileString '\n']);
                    fprintf(['Tile Count:' tileNumber '\n']);
                end
                tileCount = tileCount + 1;
                metaTileCount = metaTileCount + 1;
                
                
                tic;[qrImage] = generateQRtile(str,blockSZ);toc
                machineState = updateStateImage(machineState,r,curTile,blockSZ,rgbStates(curTile,:),satLevel(1),qrImage);
                curTile = curTile + 1;
               
                
                
                
            else
                fprintf(['****************************\n']);
                fprintf(['ERROR: Wrong tile!\n'])
                fprintf(['Please scan again.\n'])
                fprintf(['****************************\n']);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%
        % scan until the end code
        %%%%%%%%%%%%%%%%%%%%%%%%
        while ~containKVP(str,{'{command_sessionEnd}'})
            str = getValidRead(endPrompt);
            if ~containKVP(str,{'{command_sessionEnd}'})
                fprintf(['****************************\n']);
                fprintf(['Please scan the last tile!\n']);
                fprintf(['****************************\n']);
            else
                tic;[qrImage] = generateQRtile(str,blockSZ);toc
                machineState = updateStateImage(machineState,r,curTile,blockSZ,rgbStates(curTile,:),satLevel(1),qrImage);
                curTile = curTile + 1;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%
        

        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % attach the sample specific data
        fileString = [fileString sampleFileData];
        CMD = ['captureSingleImage ' deviceType ' ' usbPort ' ' ...
                folderString ' ' fileString ' ' resolution];

        fprintf(['****************************\n']);
        fprintf(['Starting data acquisition\n']);
        fprintf(['****************************\n']);
        fprintf(['Folder:' folderString '\n']); 
        fprintf(['File:' fileString '\n']);
        fprintf(['****************************\n']);
        CMD
        
        

    catch ME
        getReport(ME)
    end
end



