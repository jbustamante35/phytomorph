function [msg,qrCropBox] = getQRcode_2(I,reSZ,saveFlag)
    try
        if isdeployed
            fprintf(['*******************************\n']);
            pwd
            fprintf(['*******************************\n']);
            javaaddpath([pwd filesep 'core-3.2.1.jar']);
            javaaddpath([pwd filesep 'javase-3.2.1.jar']);
            javaclasspath
            fprintf(['*******************************\n']);
            ls
            fprintf(['*******************************\n']);
        end
        
        
        
        
        defaultAreaMinSize = 2000;
        internalBoxExpansionSize = 50;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % read image if I is char
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        saveName = '';
        if ischar(I);saveName = I;I = imread(I);end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set the resize feature to 1 as default
        if nargin == 1
            reSZ = 1;
        end
        % set the msg to nothing
        msg = '';
        % set the orginal image
        oI = I;
        % resize if requested
        if reSZ ~= 1
            I = imresize(I,reSZ);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% start QR data gather
        tm = clock;
        Lab = rgb2lab(I);
        fprintf(['Converting RGB to LAB done in:' num2str(etime(clock,tm)) '\n']);
        M = Lab(:,:,2) > 25.5;
        M = imclose(M,strel('square',51));
        M = bwlarge(imfill(M,'holes')) & M;
        M = bwlarge(M);
        R = regionprops(M);
        box = R(1).BoundingBox;
        pad = round(50*reSZ);
        box(1:2) = box(1:2) - pad;
        box(3:4) = box(3:4) + 2*pad;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find the boxes within the QR sheet
        fM = imfill(M,'holes');
        % remove small objects
        fM = bwareaopen(fM,defaultAreaMinSize*reSZ);
        % expand the internal boxes by 50 pixels
        fM = imdilate(fM,ones(round(internalBoxExpansionSize*reSZ)));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find the largest object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        R = regionprops(fM,'Area','BoundingBox');
        [J,midx] = max([R.Area]);
        if isempty(R(midx).BoundingBox)
            return
        end
        qrCropBox = R(midx).BoundingBox;

        
        
        % crop out the dataCube from the image
        dataCube = imcrop(I,qrCropBox);
        % crop out the LAB version of the data cube
        dataCubeLAB = rgb2lab(dataCube);
        %dataCubeLAB = imcrop(Lab,qrCropBox);
        % crop out the dataCube from the original size image
        dataORG = imcrop(oI,reSZ^-1*qrCropBox);
        % crop out from the masked image
        dataCubeMask = imcrop(M,qrCropBox);

        
        
        
        
        
        % added to fill in the gaps from the bad threshold values
        % note: threshold values may not close the gaps on the red frame
        % this was added at the phenome 2018 conference
        dataCubeMask = bwareaopen(dataCubeMask,round(200*reSZ));
        % added to handle text in red box - after phenome 18
        dataCubeMask = imclose(dataCubeMask,strel('disk',round(31*reSZ),0));


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find the third largest object
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dataCubeMaskf = imfill(dataCubeMask,'holes') - dataCubeMask;
        Rc = regionprops(logical(dataCubeMaskf),'Area','BoundingBox');
        [J,midxC] = sort([Rc.Area],'descend');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the qr data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        qrCube = imcrop(dataORG,(reSZ^-1)*Rc(midxC(3)).BoundingBox);
        qrCube = rgb2gray(qrCube);


        % the following set of commands are not robust
        % but found to help in X number of cases
        qrCube = imresize(qrCube,2);
        qrCube = double(localcontrast(single(qrCube)));
        if max(qrCube) > 1;qrCube = qrCube / 255;end
        qrCube = imadjust(qrCube);
        qrCube = imsharpen(qrCube,'Amount',2);
        rotValue = linspace(0,360,360);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % try to read the QR at one degree increments
        % break if read
        for e = 1:numel(rotValue)
            qrCube_read = imrotate(qrCube,rotValue(e));
            try
                msg = decode_qr(qrCube_read);
            catch ME
                fprintf(['*******************************\n']);
                getReport(ME)
                fprintf(['*******************************\n']);
                msg = [];
            end
            if ~isempty(msg)
                break
            end
        end
        % replace the newline, hardline, and tab etc
        msg = strrep(msg,char(10),'');
        msg = strrep(msg,char(13),'');
        msg = strrep(msg,char(9),'');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the day cube image
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dayCube = imcrop(dataCube,Rc(midxC(1)).BoundingBox);
        dayCubeLAB = imcrop(dataCubeLAB,Rc(midxC(1)).BoundingBox);
        %dayCubeLAB(:,:,3) = bindVec(dayCubeLAB(:,:,3));
        %dayCubeLAB(:,:,3) = adapthisteq(dayCubeLAB(:,:,3));
        dayCubeLAB(:,:,3) = dayCubeLAB(:,:,3).*dayCubeLAB(:,:,1).^-1;
        dayCubeLAB(:,:,3) = dayCubeLAB(:,:,1);
        %dayCubeLAB(:,:,3) = bindVec(dayCubeLAB(:,:,3));
        tmp = dayCubeLAB(:,:,3);
        %gm = fitgmdist(tmp(:),2);
        %[~,sIDX] = min(gm.mu);
        %BLUE = reshape(gm.cluster(tmp(:)),size(tmp)) == sIDX;
        BLUE = dayCubeLAB(:,:,3) < graythresh(dayCubeLAB(:,:,3)/100)*100;
        BLUE = bwareaopen(BLUE,round(1000*reSZ^2));
        BLUE = bwareaopen(imfill(BLUE,'holes'),round(3000*reSZ^2)) & BLUE;
        %BLUE = imclose(BLUE,strel('square',round(8*reSZ)));
        BLUE = imclearborder(BLUE);
        BLUE = bwlarge(BLUE,24);
        squareMask = BLUE;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % the below may-be good code for isolating features in the dataCube
        % on Jan, 02 2020 I am choosing to leave it "inplace"
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %imshow(squareMask,[]);
        %{
        V = rgb2hsv_fast(dayCube,'single','V');
        H = rgb2hsv_fast(dayCube,'single','H');
        % change to .50 from .6 for carrot
        squareMask = H > .45 & V < .78;
        squareMask = bwareaopen(squareMask,100);
        squareMask = imclose(squareMask,strel('disk',5));
        squareMask = imfill(squareMask,'holes');
        BOXMask = edge(squareMask);
        squareMask = imclearborder(squareMask);
        blueMask = H > .53 & H < .91;
        blueMask = bwareaopen(blueMask,400);
        blueMask = imclose(blueMask,strel('disk',5));
        %}

        %{

        %blueMask = imdilate(blueMask,strel('disk',3));

        blueMaskS = bwmorph(blueMask,'skeleton',inf);
        blueMaskS = imdilate(blueMaskS,strel('disk',2));
        % find vertical linesV
        [H, theta, rho] = hough(blueMaskS,'Theta',linspace(-10,10,20));
        P  = houghpeaks(H,12,'Threshold',0);
        linesV = houghlines(blueMaskS,theta,rho,P,'FillGap',200,'MinLength',250);
        % find horizontal linesH
        [H, theta, rho] = hough(blueMaskS','Theta',linspace(-10,10,20));
        P  = houghpeaks(H,8,'Threshold',0);
        linesH = houghlines(blueMaskS',theta,rho,P,'FillGap',200,'MinLength',600);

        %{
        imshow(dayCube,[]);
        hold on
        for k = 1:length(linesV)
           xy = [linesV(k).point1; linesV(k).point2];
           plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
           % Plot beginnings and ends of linesV
           plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
           plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
        end
        for k = 1:length(linesH)
           xy = [fliplr(linesH(k).point1); fliplr(linesH(k).point2)];
           plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
           % Plot beginnings and ends of linesV
           plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
           plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
        end
        %}

        BOXMask = zeros(size(blueMask));
        for k = 1:length(linesV)
           xy = [linesV(k).point1; linesV(k).point2];
           X = round(linspace(xy(1,1),xy(2,1),3000));
           Y = round(linspace(xy(1,2),xy(2,2),3000));
           for l = 1:numel(X)
               BOXMask(Y(l),X(l)) = 1;
           end
        end
        for k = 1:length(linesH)
           xy = [fliplr(linesH(k).point1); fliplr(linesH(k).point2)];
           X = round(linspace(xy(1,1),xy(2,1),3000));
           Y = round(linspace(xy(1,2),xy(2,2),3000));
           for l = 1:numel(X)
               BOXMask(Y(l),X(l)) = 1;
           end
        end

        %}


        %{
        BOXMask = imdilate(BOXMask,strel('disk',7));
        blueMask = blueMask.*BOXMask;

        blueMask = imdilate(blueMask,strel('disk',3));
        markerStrip = V < .6 - blueMask;
        markerStrip = bwareaopen(markerStrip,100);
        %}

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the centers of each day cube
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ERO = round(50*reSZ);
        S1 = std(imfill(squareMask,'holes'),1,1);
        S1 = S1 + .0001*rand(size(S1));
        %bk1 = imerode(S1,ones(1,ERO));
        %bk1 = imfilter(bk1,fspecial('disk',51),'replicate');
        S1 = imfilter(S1,fspecial('disk',round(31*reSZ)),'replicate');
        [lm1] = nonmaxsuppts(S1,round(101*reSZ));
        R1 = regionprops(lm1,'Centroid');
        lm1 = zeros(size(lm1));
        for e = 1:numel(R1)
            lm1(round(R1(e).Centroid(1))) = 1;
        end
        S1 = bindVec(S1);
        level = graythresh(S1);
        lm1 = lm1 & S1 > level;

        S2 = std(imfill(squareMask,'holes'),1,2);
        S1 = S2 + .0001*rand(size(S1));
        %bk2 = imerode(S2,ones(ERO,1));
        %bk2 = imfilter(bk2,fspecial('disk',51),'replicate');
        S2 = imfilter(S2,fspecial('disk',round(31*reSZ)),'replicate');
        [lm2] = nonmaxsuppts(S2,round(81*reSZ));
        R2 = regionprops(lm2,'Centroid');
        lm2 = zeros(size(lm2));
        for e = 1:numel(R2)
            lm2(round(R2(e).Centroid(2))) = 1;
        end
        S2 = bindVec(S2);
        level = graythresh(S2);
        lm2 = lm2 & S2 > level;

        M = double(lm2)*double(lm1);
        [c1 c2] = find(M);
        centers = [c1 c2];


        %{
        imshow(dayCube,[]);
        hold on
        plot(c2,c1,'*');
        for e = 1:numel(c1)
            plot(c2(e),c1(e),'*');
            text(c2(e),c1(e),num2str(e));
        end
        %}


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % snap each center to nearest centroid 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        R = regionprops(squareMask,'Centroid','BoundingBox');
        CEN = [];
        for e = 1:numel(R)
            CEN = [CEN;fliplr(R(e).Centroid)];
        end
        for e = 1:size(centers,1)
            idx = snapTo(CEN,centers(e,:));
            boxTmp = imcrop(squareMask,R(idx).BoundingBox);
            QUARTER = R(idx).BoundingBox(3:4)/4;
            HALF = R(idx).BoundingBox(3:4)/2;
            newBOX = R(idx).BoundingBox;
            newBOX = newBOX(1:2) + QUARTER;
            newBOX(3:4) = HALF;
            boxCenter = imcrop(squareMask,newBOX);
            %markerTmp = imcrop(markerStrip,R(idx).BoundingBox);
            %checked(e) = sum(boxTmp(:)-markerTmp(:))/sum(boxTmp(:)) < .98;
            checked(e) = sum(boxCenter(:)) > 100*reSZ^2;
        end
        dayMatrix = [4:27];
        didx = find(checked);
        if ~isempty(didx)
            pictureDay = dayMatrix(didx(end));
        else
            pictureDay = 0;
        end
        msg = [';' msg ';PictureDay:' num2str(pictureDay) ';'];

        msg = strrep(msg,';',';;');
        fidx = strfind(msg,';');
        msg(fidx(1:2:end)) = '}';
        msg(fidx(2:2:end)) = '{';
        msg(1) = [];
        msg(end) = [];
        msg = strrep(msg,':','_');
        msg = strrep(msg,'/','-');
        
        
        
        if ~isempty(saveName)
            fprintf(['**********************\n']);
            fprintf(['Entering saving phase:\n']);
            fprintf(['**********************\n']);
            fprintf(['MSG:' msg '\n']);
            fprintf(['**********************\n']);
            
            [~,saveName] = fileparts(saveName);
            oPath = './output/';
            mkdir(oPath);
            oFile = [oPath saveName];
            fileID = fopen(oFile,'w');
            fprintf(fileID,'%s',msg);
            fclose(fileID);
        end
    catch ME
        msg = '';
        qrCropBox = [];
    end
end

%{


    close all
    FilePath = '/mnt/tetra/nate/projectData/maizeSeedling/NEF_to_tiff/';
    FileList = {};
    FileExt = {'tif'};
    FileList = gdig(FilePath,FileList,FileExt,1);

    bad = [];
   
 



    bad = [];
    bad2 = [];
    bCheck = zeros(numel(FileList),1);
    parfor e = 1:numel(FileList)

      
            e
            tic
            I = imread(FileList{e});
            % strip off the right side of the image
            I(:,end-60:end,:) = [];
            % get the file parts
            [pth,nm,ext] = fileparts(FileList{e});
            [msg{e},qrCropBox{e}] = getQRcode_2(I,.35);

            bCheck(e) = strcmp(nm,msg{e});


            toc

    end


    %{
    CL = reshape(I,[size(I,1)*size(I,2) 3]);
    hardMean = [0.8249    0.4369    0.3684];
    hardCov = [[0.0085    0.0094    0.0087];[0.0094    0.0109    0.0101];[0.0087    0.0101    0.0095]];
    P = mvnpdf(CL,hardMean,hardCov);
    P = reshape(P,[size(I,1) size(I,2)]);
    M = log(P) > -20;
    %}
    %{
    H = rgb2hsv_fast(I,'single','H');
    hsvI = rgb2hsv(I);
    H = hsvI(:,:,1);
    % also changed the threshold value - OLD .04
    val = .044;
    M = H < val | H > (1 - val);
    %}

%}