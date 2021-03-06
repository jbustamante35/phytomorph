function [finalScore,nX,rawX] = detectEmergence_ver2(FileList,emergenceNet,frameNet,nsz,zU,zE,CELLMASK,BOX,vidx,Tscore,GMM)
    try

        % start clock
        tm = clock;
        mkdir('./output/')
        fprintf(['Reporting on GMM:\n']);
        GMM
        fprintf(['Reporting on GMM:\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% STARTING EMERGENCE NET
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isdeployed()
            poolObj = parpool(4);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % remove non-tifs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isdeployed()
            fprintf(['******************************************************\n']);
            fprintf(['starting: removing non tiffs: \n']);tic
            fprintf(['******************************************************\n']);
            for e = 1:numel(FileList)
                [pth,nm,ext] = fileparts(FileList{e});
                if ~strcmp(lower(ext(1:4)),'.tif')
                    rm(e) = true;
                else
                    rm(e) = false;
                end
            end
            FileList(find(rm)) = [];
            fprintf(['******************************************************\n']);
            fprintf(['ending: removing non tiffs:' num2str(toc) '\n']);
            fprintf(['******************************************************\n']);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  remove non-tifs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sort the image list
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['******************************************************\n']);
        fprintf(['starting: file sort: \n']);tic
        fprintf(['******************************************************\n']);
        for e = 1:numel(FileList)
            [pth,nm,ext] = fileparts(FileList{e});
            NAME(e) = str2num(nm);
        end
        [~,sidx] = sort(NAME);
        FileList = FileList(sidx);
        framesToMeasure = min(numel(FileList),250);
        framesToMeasure = numel(FileList);
        FileList = FileList(1:framesToMeasure);
        fprintf(['******************************************************\n']);
        fprintf(['ending: file sort:' num2str(toc) '\n']);
        fprintf(['******************************************************\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sort the image list
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % replace bad images
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isdeployed()
            fprintf(['******************************************************\n']);
            fprintf(['starting: removing corrupt images: \n']);tic
            fprintf(['******************************************************\n']);
            for e = 1:numel(FileList)
                try
                    tmpI = imread(FileList{e});
                    good(e) = true;
                catch
                    good(e) = false;
                end
            end
            IDX = 1:numel(FileList);
            bidx = find(good==false);
            gidx = find(good==true);
            fprintf(['Found good:' num2str(numel(gidx)) ':Found bad:' num2str(numel(bidx)) '\n']);
            for e = 1:numel(bidx)
                [J,ridx] = min(abs(gidx - bidx(e)));
                IDX(bidx(e)) = gidx(ridx);
                fprintf(['corrupt image @:' num2str(bidx(e)) '<-->' num2str(gidx(ridx)) '\n']);
            end
            FileList = FileList(IDX);
            fprintf(['******************************************************\n']);
            fprintf(['ending: removing corrupt images:' num2str(toc) '\n']);
            fprintf(['******************************************************\n']);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % replace bad images
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % reporting
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['Measuring total number of frames:' num2str(numel(FileList)) '\n']);
        fprintf(['******************************************************\n']);
        fprintf('Reporting first two and last two file names\n');
        FileList{1}
        FileList{2}
        FileList{end-1}
        FileList{end}
        fprintf('Reporting first two and last two file names\n');
        fprintf(['******************************************************\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % reporting
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get camera parameters from the checkerboard
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['******************************************************\n']);
        fprintf(['starting: get checker board \n']);tic
        fprintf(['******************************************************\n']);
        [tform] = getCameraPara(FileList{1});
        fprintf(['Reporting t_form:\n']);
        tform.T
        fprintf(['Reporting t_form:\n']);
        fprintf(['******************************************************\n']);
        fprintf(['ending: get checker board in:' num2str(toc) '\n']);
        fprintf(['******************************************************\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get camera parameters from the checkerboard 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        %[boundingBox,centerPoints,MASK,I,uCL,covCL] = getCropBoxes(FileList,tform,168,10,150,[],false);  
      


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % register the images
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        toRegister = false;
        if toRegister
            fprintf(['******************************************************\n']);
            fprintf(['starting: image registration \n']);tic
            fprintf(['******************************************************\n']);
            [cameraShift] = registerAgainstBackground(FileList,800,0,tform,uCL,covCL);
            fprintf(['******************************************************\n']);
            fprintf(['ending: image registration:' num2str(toc) '\n']);
            fprintf(['******************************************************\n']);
        else
            cameraShift = [];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  register the images
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get crop boxes and order them
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['******************************************************\n']);
        fprintf(['starting: get crop boxes \n']);tic
        fprintf(['******************************************************\n']);
        [boundingBox,centerPoints,MASK,I] = getCropBoxes(FileList,tform,168,10,150,cameraShift,false,GMM);
        fprintf(['******************************************************\n']);
        fprintf(['ending: get crop boxes:' num2str(toc) '\n']);
        fprintf(['******************************************************\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['******************************************************\n']);
        fprintf(['starting: order crop boxes \n']);tic
        fprintf(['******************************************************\n']);
        [boundingBox,centerPoints,LABELS] = orderCropBoxes(MASK,boundingBox,centerPoints,I,false);
        fprintf(['******************************************************\n']);
        fprintf(['ending: order crop boxes:' num2str(toc) '\n']);
        fprintf(['******************************************************\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get crop boxes and order them
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %tmp = getRectifiedImage(FileList{229},tform,cameraShift{229});
        %J = imcrop(tmp,boundingBox{155});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % crop and measure each cell
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [finalScore,nX,rawX] = cropTimeSeriesFromRectified_ver3(FileList,boundingBox,[1:1:framesToMeasure],tform,'',emergenceNet,frameNet,nsz,zU,zE,CELLMASK,BOX,vidx);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % crop and measure each cell
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % disk spool
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isdeployed()
            fprintf(['Starting disk spool.\n']);
            mkdir('./output/')
            [pth,nm,ext] = fileparts(FileList{1});
            fidx = strfind(pth,filesep);
            outName = pth;
            csvwrite(['./output/' outName '.csv'],finalScore');
            save(['./output/' outName '_dataPackage.mat'],'finalScore','nX','rawX');
            
            for n = 1:numel(FileList)
                [pth,nm,ext] = fileparts(FileList{n});
                numberName(n) = str2num(nm);
            end
            
            for n = 1:numel(finalScore)
                if finalScore(n) ~= 0
                    finalScore_frame(n) = numberName(finalScore(n));
                else
                    finalScore_frame(n) = 0;
                end
            end
            csvwrite(['./output/' outName '_fileName.csv'],finalScore_frame');
            
            fprintf(['Ending disk spool.\n']);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % disk spool
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        

    catch ME
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % clean up
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isdeployed()
            %fprintf('Shutting down par pool.\n');
            %delete(poolObj);
            %{
            poolObj = gcp('nocreate');
            while ~isempty(poolObj)
                delete(poolObj);
                poolObj = gcp('nocreate');
            end
            %}
        end
        fprintf(['Starting close of all figures.\n']);
        close all
        fprintf(['Ending close of all figures.\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % clean up
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        getReport(ME)
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % clean up
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isdeployed()
        %fprintf('Shutting down par pool.\n');
        %poolObj = gcp('nocreate')
        %while ~isempty(poolObj)
        %    delete(poolObj)
        %    poolObj = gcp('nocreate')
        %end
        %clear poolObj
    end
    fprintf(['Starting close of all figures.\n']);
    close all
    fprintf(['Ending close of all figures.\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % clean up
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % report total run time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['Done running TOTAL time:' num2str(etime(clock,tm)) '\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % report total run time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
        
end


%{
    FilePath = '/mnt/tetra/nate/cornPopper/';
    FileExt = {'tiff'};
    FileList = gdig(FilePath,{},FileExt,1);
    obj.func(FileList(1:200))


    FilePath = '/mnt/tetra/nate/junk/20180601_Rack1_Camera1/';
    FileExt = {'tiff'};
    FileList = gdig(FilePath,{},FileExt,1);
    obj.func(FileList)
    

%}