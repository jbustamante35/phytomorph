function [] = deviceBank(debug)
    try
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
        fprintf(['Session started for ' deviceType '@' usbPort ' usbPort.\n']);
        expectedTilesN = totalCodeN;
        %%%%%%%%%%%%%%%%%%%%%%%%

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
                sampleCodeN = sampleCodeN - 1;
            else
                fprintf('Your command was not a valid sample QR.\n');
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%
        % scan meta data codes until the top is scanned
        %%%%%%%%%%%%%%%%%%%%%%%%
        str = '';
        flag=false;
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
            else
                fprintf(['****************************\n']);
                fprintf(['ERROR: Wrong tile!\n'])
                fprintf(['Please scan again.\n'])
                fprintf(['****************************\n']);
            end
        end



        %%%%%%%%%%%%%%%%%%%%%%%%
        % scan until the end code 
        str = '';
        while ~containKVP(str,{'{command_sessionEnd}'})
            str = getValidRead(endPrompt);
            if ~containKVP(str,{'{command_sessionEnd}'})
                fprintf(['****************************\n']);
                fprintf(['Please scan the last tile!\n']);
                fprintf(['****************************\n']);

            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%
        % attach the sample specific data
        fileString = [fileString sampleFileData];


        fprintf(['****************************\n']);
        fprintf(['Starting data aquizition\n']);
        fprintf(['****************************\n']);
        fprintf(['Folder:' folderString '\n']); 
        fprintf(['File:' fileString '\n']);
        fprintf(['****************************\n']);

    catch ME
        ME;
    end
end












