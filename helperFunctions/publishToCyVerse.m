function [] = publishToCyVerse(pubtype)
    if nargin == 0
        
        functionFile = '-a smartMain_new_ver1.m -a singleGreenFieldImage.m -a inpoly.m -a heavisideTransitionFunction.m -a folderTIPS.m -a boo.mexa64 -a decode_qr.m -a mj_getFields.m -a getPlantCell_and_Mask.m -a mergeCropBoxes.m -a getPlantBoxesForFile.m -a checkerBoardAnalysis.m -a getFieldForFileName.m -a getCheckBoardfileNames.m -a getDayfileNames.m -a parseAndLabelStack.m -a overHead_main.m -a generalizeFeatureExtractor.m -a generalizeLoader.m -a applyAllLayers.m -a trainAIlayer.m -a cRunner.m -a scanalyzerMain.m -a /usr/local/MATLAB/R2017b/toolbox/images/imdata/eSFRdefaultGrayReference.mat -a /usr/local/MATLAB/R2017b/toolbox/images/imdata/eSFRdefaultColorReference.mat -a cellForDirk.m -a generatePlant.m -a generateMetaDataTree.m -a stomataMaster.m -a fftPatch.m -a myInterp2Sampler.m -a applyNetworkToPoint.m -a generateImageDomain.m -a detectEmergence_ver2.m -a rgbAndhsvExtract.m -a arborist.m -a measureCrossOver.m -a extractPhenotypesFromOverheadCamera_ver2.m -a clusterForest.m -a clusterTree.m -a clusterNode.m -a extractPhenotypesFromOverheadCamera.m -a genQRLargeFormatSheets.m -a detectEmergence -a imsubtract.m -a rgb2hsv_fast.m -a generateImageClass.m -a getPNN_func.m -a shapeVerticalStripNozzle.m -a network.m -a func_depthStack.m -a func_resizeDepthStack.m -a constantTransitionFunction.m -a myProb.m -a my_hmm.m -a hmm_node.m -a nozzleManifold.m -a cropImages_v2.m -a func_thumbNail.m -a smartMain_v4.m -a smartMain_v2.m -a sorguhmLeafAnalysis.m -a maizeSeedling_func3.m -a maizeSeedling_func2.m -a maizeSeedling_func1.m -a smartMain.m -a singleWholeCarrotStage2.m -a petLength.m -a partialFunction.m -a condorDeploy_ver0.m -a confocalRatio.m -a isolateRoots_overStack.m -a mecka.m -a singleWholeCarrotAnalyze.m -a op0.m -a singleSeedlingImage.m';
        cdir = dir('/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/Maize/maizeSeedling/nets/');
        cdir(1:2) = [];
        for e = 1:numel(cdir)
            functionFile = ['-a ' cdir(e).name ' ' functionFile];
        end
        tmpCompileDirectory = '/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/deDeployR2017a/';
        mkdir(tmpCompileDirectory)
        tmpCompileDirectory = '/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/deDeploy/';
        mkdir(tmpCompileDirectory)

        CMD = ['mcc -d ' tmpCompileDirectory ' ' functionFile ' -a  gmdistribution.m -a cJob.m -a im2single.m -m -v DEwrapper.m'];
        eval(CMD);
        pushCMD = ['iput -f /mnt/scratch1/phytomorph_dev/Extraction/Pipelines/deDeploy/DEwrapper /iplant/home/nmiller/publicData/DEwrapper'];
        [pushR] = system(pushCMD,'-echo');
        
        CMD = ['mcc -d ' tmpCompileDirectory ' ' functionFile ' -a  gmdistribution.m -a cJob.m -a im2single.m -m -v DEwrapper_osg.m'];
        eval(CMD);
        pushCMD = ['iput -f /mnt/scratch1/phytomorph_dev/Extraction/Pipelines/deDeploy/DEwrapper_osg /iplant/home/nmiller/phytoMorph_public_deploy/matlab/DEwrapper'];
        [pushR] = system(pushCMD,'-echo');
        
        
        
        
        
    elseif nargin == 1
        switch pubtype
            case 'algorithmRouterTable'
                %CMD = ['iput -f /mnt/snapper/nate/phytoKeyPool/algoRouteTable.mat /iplant/home/nmiller/publicData/']
                CMD = ['iput -f  /mnt/snapper/nate/phytoKeyPool/algoRouteTable.json /iplant/home/nmiller/publicData/'];
                [~,r] = system(CMD);
                r
            
            case 'plantProfiler'
                D = datestr(datetime,'ss_hh_dd_mm_yyyy');
                 
                basePath = '/mnt/scratch1/phytomorph_dev/Deploy/';
                iBasePath = '/iplant/home/nmiller/publicData/';
                iWholePath = [iBasePath pubtype filesep];
                CMD = ['imkdir -p ' iWholePath];
                
                wholePath = [basePath pubtype filesep];
                mmkdir(wholePath);
                
                
                programFile = [wholePath pubtype];
                runScript = [wholePath 'run_' pubtype '.sh'];
                
                iprogramFile = [iWholePath pubtype];
                irunScript = [iWholePath 'run_' pubtype '.sh'];
                
                mcc -m /mnt/scratch1/phytomorph_dev/Aquire/plantProfiler/plantProfiler.m ...
                      -d /mnt/scratch1/phytomorph_dev/Deploy/plantProfiler
                fprintf(['Backing up old copy.\n']);
                CMD = ['icp -f ' iprogramFile ' ' iprogramFile D];
                system(CMD);
                fprintf(['Updating new copy.\n']);
                CMD = ['iput -f ' programFile ' ' iprogramFile];
                system(CMD);
                fprintf(['Updating shell run script.\n']);
                CMD = ['iput -f ' runScript ' ' irunScript];
                system(CMD);
                
            case 'deviceBank'
                D = datestr(datetime,'ss_hh_dd_mm_yyyy');
                 
                basePath = '/mnt/scratch1/phytomorph_dev/Deploy/';
                iBasePath = '/iplant/home/nmiller/publicData/';
                iWholePath = [iBasePath pubtype filesep];
                CMD = ['imkdir -p ' iWholePath];
                
                wholePath = [basePath pubtype filesep];
                mmkdir(wholePath);
                
                
                programFile = [wholePath pubtype];
                runScript = [wholePath 'run_' pubtype '.sh'];
                
                iprogramFile = [iWholePath pubtype];
                irunScript = [iWholePath 'run_' pubtype '.sh'];
                
                localD = '/mnt/scratch1/phytomorph_dev/Deploy/deviceBank';
                mkdir(localD);
                
                % pay attention here to the strings with copying
                mcc -m /mnt/scratch1/phytomorph_dev/Aquire/deviceBank/deviceBank.m ...
                      -d /mnt/scratch1/phytomorph_dev/Deploy/deviceBank
                
                fprintf(['Backing up old copy.\n']);
                CMD = ['icp -f ' iprogramFile ' ' iprogramFile D];
                system(CMD);
                fprintf(['Updating new copy.\n']);
                CMD = ['iput -f ' programFile ' ' iprogramFile];
                system(CMD);
                fprintf(['Updating shell run script.\n']);
                CMD = ['iput -f ' runScript ' ' irunScript];
                system(CMD);
            
            case 'phytoJSONdictionary'
                pmd = generatePhytomorphDictionary();
                pmd = jsonencode(pmd);
                localCopyFile = '~/phytoMorphTK/baseInstall/pmd.json';
                [pth,nm,ext] = fileparts(localCopyFile);
                mmkdir(pth);
                % save local
                fileID = fopen(localCopyFile,'w');
                fprintf(fileID,'%s',pmd);
            case 'gLoader'
                %%%%%%%%%%%
                % to publish generalized loader
                loaderFunc = @(X,loaderType,loaderArgs)generalizeLoader(X,loaderType,loaderArgs);
                loaderWrapper = partialFunction(loaderFunc,'generalizeLoader');
                loaderWrapper.publish();
             case 'gFeatureExtractor'
                %%%%%%%%%%%
                % to publish generalized loader
                extractFunction = @(X,extractType,extractArgs)generalizeFeatureExtractor(X,extractType,extractArgs);
                extractFunction = partialFunction(extractFunction,'generalizeFeatureExtractor');
                extractFunction.publish();
        end
    end
end





