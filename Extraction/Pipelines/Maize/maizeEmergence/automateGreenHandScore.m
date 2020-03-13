function [greenPoppedFrame] = automateGreenHandScore(useFlag,pathList,imageScale2,centerMask,storePath)


    toAVG = 50;
    dropThresh = .5;
    greenThresh = 1;
    greenPoppedFrame = [];
    startTM = clock;
    parfor gidx = 1:numel(pathList)
        if useFlag(gidx)
            try
                tic;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                % scan and load the data
                mFilePath = pathList{gidx};
                mFileList = {};
                FileExt = {'tif'};
                mFileList = sdig(mFilePath,mFileList,FileExt,1);
                mFileList = mFileList{1};
                tkp = [];
                for e = 1:numel(mFileList)
                    tkp(e) = contains(mFileList{e},'raw');
                end
                mFileList = mFileList(logical(tkp));
                NM = [];
                for e = 1:numel(mFileList)
                    [pth,nm,ext] = fileparts(mFileList{e});
                    fidx = strfind(nm,'_');
                    nm = nm(1:(fidx(1)-1));
                    NM(e) = str2num(nm);
                end
                [~,sidx] = sort(NM);
                mFileList = mFileList(sidx);

                tLABsig = [];
                %tLABsig2 = [];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                for e = 1:numel(mFileList)
                    tmpI = double(imread(mFileList{e}))/255;
                    tmpI = imresize(tmpI,imageScale2);
                    tmpI = bsxfun(@times,tmpI,centerMask);
                    tmpLab = rgb2lab(tmpI);
                    %imshow(tmpLab,[]);
                    %drawnow
                    szLAB = size(tmpLab);
                    tLABsig(e,:) = mean(reshape(tmpLab,[prod(szLAB(1:2)) szLAB(3)]),1);
                    %tLABsig2(e,:) = std(reshape(tmpLab,[prod(szLAB(1:2)) szLAB(3)]),1,1);
                end



                sig = tLABsig(:,2);
                sig = imfilter(sig,fspecial('average',[11 1]),'replicate');
                baseLine = mean(sig(1:toAVG));
                sig = sig - baseLine;
                greenPulse = sig < -dropThresh;
                greenDELTA = min(sig);
                isGreen = greenDELTA < -greenThresh;

                gridx = find(greenPulse);

                if ~isempty(gridx) & isGreen
                    greenPopped = gridx(1);
                else
                    greenPopped = 0;
                end

                greenPoppedFrame(gidx) = greenPopped;


                perTM = toc;
                TOTALesti = perTM*numel(pathList);
                goneTM = etime(clock,startTM);
                leftTM = (TOTALesti - goneTM)/60/12;


                fprintf(['Done with automated-hand green measure:' num2str(gidx) ':' num2str(numel(pathList))  ...
                    ':' num2str(perTM/60) 'min :' num2str(leftTM)  'min. \n']);
            catch ME
                ME
            end
        else
            greenPoppedFrame(gidx) = -1;
        end
    end
    
    
    
    FilePath = storePath;
    csvFileList = {};
    FileExt = {'csv'};
    csvFileList = fdig(FilePath,csvFileList,FileExt,1);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % scan for only csv files that have the key word in the name
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kp = [];
    for e = 1:numel(csvFileList)
        [p,n] = fileparts(csvFileList{e});
        kp(e) = contains(n,'Handscore') & ~contains(n,'BK') & ~contains(n,'BAD');
    end
    csvFileList = csvFileList(find(kp));
    csvFile = csvFileList{1};
    csvFile = strrep(csvFile,'Handscore','Leafscore');
    
    csvwrite(csvFile,greenPoppedFrame');



end