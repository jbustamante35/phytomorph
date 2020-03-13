classdef networkTester
    
    
    properties (Static)
        testProcedureRoot = '/mnt/scratch1/myLinks/octerine/testProcedures/';
        testProcedureData = 'data/';
        testProecdureTest = 'tests/';
    end
    
    
    properties
        
        dut; % device under test
        
    end
    
    methods
        
        function [] = setDUT(obj,dut)
            obj.dut = dut;
        end
        
        
    end
    
    methods (Static)
        

 


        function [p] = dataPath()
            p = [networkTester.testProcedureRoot networkTester.testProcedureData];
        end
        
        function [p] = testPath()
            p = [networkTester.testProcedureRoot networkTester.testProecdureTest];
        end
        
        function [p] = procedureTestPath(pType,spType,testDate,testVersion)
            p = [networkTester.testPath pType filesep spType filesep '{testVersion_' testVersion '}' '{testDate_' testDate '}' filesep];
        end

        function [newFileList] = addImageToTestProcedure(sourceList,network,spType,sType)
            pType = network.testName;
            % pType = procedure type
            % spType = sub procedure type(s)
            % sType = link or copy
            if ~iscell(spType);spType = {spType};end
            for s = 1:numel(sourceList)
                target = {};
                source = sourceList{s};
                [pth,nm,ext] = fileparts(source);
                for e = 1:numel(spType)
                    sourceDirectory = [networkTester.dataPath ...
                        pType filesep spType{e} filesep];
                    CMD = ['mkdir -p "' sourceDirectory '"'];
                    [r,o] = system(CMD);
                    target{e} =  [sourceDirectory hash(source) '_' nm ext];
                end


                for e = 1:numel(target)
                    fprintf([source '-->' target{e} '\n']);
                    switch sType
                        case 'copy'
                            if ~exist(target{e},'file')
                                copyfile(source,target{e});
                            end
                        case 'link'
                            CMD = ['ln -s ''' source ''' ''' target{e} ''''];
                            [r,o] = system(CMD);
                    end
                end
                newFileList{s} = target;
            end
          
        end
        
        function [FileList] = listImagesForTestProcedure(network,spType)
            pType = network.testName;
            FilePath = [networkTester.dataPath ...
                        pType filesep spType filesep];
            FileList = {};
            FileExt = {'TIFF','TIF','JPG','TIF','TIFF','PNG','BMP'};
            FileExt = cat(2,FileExt,lower(FileExt));
            FileList = gdig(FilePath,FileList,FileExt,1);
        end

        function [] = loopOverData(dut,FileList,testLocation)
            % for each image in the list
            for e = 1:numel(FileList)
                try
                    % clear th global timing datad
                    timingBlock('clear');

                    % extract the file from the list
                    % extract the hash from the file name
                    inputFile = FileList{e};
                    [~,nm] = fileparts(inputFile);
                    fidx = strfind(nm,'_');
                    hashV = nm(1:(fidx(1)-1));
                    oPath = [testLocation '{imageHash_' hashV '}/'];

                    % construct the input
                    inputStruct = {inputFile,oPath};


                    testFluid = digitalFluid('testFluid',inputStruct);
                    testFluid = testFluid.flow(dut);


                    if ~isempty(testFluid.errorLog)
                        here = 1;
                    end

                catch ME
                    ME;
                end


            end
        end

        function [] = test(dut,spType)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % extract the procedure name from the layername
            pType = dut.testName;
            % make test 
            testDate = datestr(datetime,'FFF_ss_hh_dd_mm_yyyy');
            % version
            standardVersion = num2str(dut.algorithmVersion);
            % build the test location
            curTestLocation = networkTester.procedureTestPath(pType,spType,testDate,standardVersion);
            % make the output location
            CMD = ['mkdir -p ''' curTestLocation ''];
            [r,o] = system(CMD);
            % list the needed data for the dut
            [FileList] = networkTester.listImagesForTestProcedure(dut,spType);
            % loop over the data
            networkTester.loopOverData(dut,FileList,curTestLocation);
            

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % compare hash of outputs
            [numCompare,diffSet] = networkTester.hashCompare(dut,spType,curTestLocation,standardVersion);
           
            standardTable = networkTester.dataCompare(dut,spType,curTestLocation,standardVersion);
        end

        function [] = setStandard(dut,spType)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % extract the procedure name from the layername
            pType = dut.testName;
            % make test 
            testDate = 'standard';
            % version
            standardVersion = num2str(dut.algorithmVersion);
            % build the test location
            curTestLocation = networkTester.procedureTestPath(pType,spType,testDate,standardVersion);
            % make the output location
            CMD = ['mkdir -p ''' curTestLocation ''];
            [r,o] = system(CMD);
            % list the needed data for the dut
            [FileList] = networkTester.listImagesForTestProcedure(dut,spType);
            % loop over the data
            networkTester.loopOverData(dut,FileList,curTestLocation);
            
            
           
        end

        function [numCompare,diffSet] = hashCompare(dut,spType,curTestLocation,standardVersion)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % scan current test directory for data files
            FilePath = curTestLocation;
            tFileList = {};
            FileExt = {'tiff','TIF','csv','mat','jpg'};
            verbose = 1;
            tFileList = fdig(FilePath,tFileList,FileExt,verbose);
            testHashTable = fileHash(tFileList,true);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % extract the procedure name from the layername
            pType = dut.testName;
            % make test 
            testDate = 'standard';
            % buld the test location
            standardLocation = networkTester.procedureTestPath(pType,spType,testDate,standardVersion);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % scan current test directory for data files
            FilePath = standardLocation;
            tFileList = {};
            FileExt = {'tiff','TIF','csv','mat','jpg'};
            verbose = 1;
            tFileList = fdig(FilePath,tFileList,FileExt,verbose);
            standardHashTable = fileHash(tFileList,true);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % sort data
            testHashTable = sortrows(testHashTable,'hashValue');
            standardHashTable = sortrows(standardHashTable,'hashValue');

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % simple file number compare
            numCompare = size(testHashTable,1) == size(standardHashTable,1);

             

            [tMissList] = ...
                compareHashedFileTables(testHashTable,standardHashTable);

            [sMissList] = ...
                compareHashedFileTables(standardHashTable,testHashTable);
            
            diffSet{1} = sMissList;
            diffSet{2} = tMissList;

        end

        function [standardTable] = dataCompare(dut,spType,curTestLocation,standardVersion)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % scan current test directory for data files
            FilePath = curTestLocation;
            cFileList = {};
            FileExt = {'csv'};
            verbose = 1;
            cFileList = sfdig(FilePath,cFileList,FileExt,verbose);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % extract the procedure name from the layername
            pType = dut.testName;
            % make test 
            testDate = 'standard';
            % buld the test location
            standardLocation = networkTester.procedureTestPath(pType,spType,testDate,standardVersion);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % scan current test directory for data files
            FilePath = standardLocation;
            tFileList = {};
            FileExt = {'csv'};
            verbose = 1;
            sFileList = sfdig(FilePath,tFileList,FileExt,verbose);
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % sort current file list
            for e = 1:numel(cFileList)
                [pth,nm,ext] = fileparts(cFileList{e}{1});
                fidx = strfind(pth,filesep);
                cToSort{e} = pth((fidx(end)+1):end);
            end
            [cToSort,sidx] = sort(cToSort);
            cFileList = cFileList(sidx);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % sort standard file list
            for e = 1:numel(sFileList)
                [pth,nm,ext] = fileparts(sFileList{e}{1});
                fidx = strfind(pth,filesep);
                sToSort{e} = pth((fidx(end)+1):end);
            end
            [sToSort,sidx] = sort(sToSort);
            sFileList = sFileList(sidx);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 
            cnt = 1;
            X = {};
            for e = 1:numel(sFileList)
                curBase = fileparts(cFileList{e}{1});
                stdBase = fileparts(sFileList{e}{1});
                for f = 1:numel(sFileList{e})


                    [pth,csvNames{cnt}] = fileparts(sFileList{e}{f});

            

                    curStdFile = sFileList{e}{f};
                    curTstFile = strrep(curStdFile,stdBase,curBase);
                    fprintf([curStdFile '<->' curTstFile '\n']);
                    standardData = cell2mat(readtext(curStdFile));
                    testData = cell2mat(readtext(curTstFile));

                    X{cnt} = [standardData(:),testData(:)];

                    if ~isempty(standardData)
                        corrValues(cnt) = corr(standardData(:),testData(:));
                        cnt = cnt + 1;
                    end
                end
            end


            for e = 1:numel(csvNames)
                try
                    [pth,nm,ext] = fileparts(csvNames{e});
                    strTag = '{fileName_';
                    strIdx = strfind(nm,strTag);
                    stpIdx = strfind(nm,'}');
                    subIdx = find(stpIdx > strIdx);
                    nm(strIdx(1):stpIdx(subIdx(1))) = [];
                    csvNames{e} = [nm];
                catch

                end
            end

            UQ = unique(csvNames);
            CV = [];
            uCV = [];
            for u = 1:numel(UQ)
                fidx = find(strcmp(csvNames,UQ{u}));
                tmpX = cell2mat(X(fidx)');
                %{
                plot(tmpX(:,1),tmpX(:,2),'.');
                drawnow
                %}
                CV(u) = corr(tmpX(:,1),tmpX(:,2));
                uCV(u) = mean(corrValues(fidx));
            end
            standardTable = table(UQ',CV',uCV','VariableNames',{'fileNames','massCorrelation','meanCorrelation'});



        end

    end


    
    
end