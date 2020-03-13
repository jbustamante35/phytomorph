
FilePath = '/home/nate/Vein_pattern_images/';
%FilePath = '/home/nate/Vein_pattern_images/OMY183_Col/';
cFileList = {};
FileExt = {'tif'};
cFileList = fdig(FilePath,cFileList,FileExt,1);
%%
for e = 1:numel(cFileList)
    I = imread(cFileList{e});
    imshow(I,[]);
    drawnow
end
%%
close all
mag = -.4;
mag = 0;






SEN = linspace(.35,.45,8);

for s = 1:numel(SEN)

    for e = 1:numel(cFileList)

        [pth,nm,ext] = fileparts(cFileList{e});
        I = double(imread(cFileList{e}))/255;
        I = rgb2gray(I);
        M = I > graythresh(I);
        M = bwareaopen(M,3000);
        M = imclearborder(M);

        R = regionprops(logical(M),'BoundingBox','Area');
        


        for b = 1:numel(R)


        % value for spuring the skeleton
        spurAmount = 3;
        % value to measure the cost for spurs that connect branch to end
        snipAmount = 11;
        % amount to dilate the skeleton
        skeletonDilateAmount = 3;
        % amount to erode the mask
        totalMaskErodeAmount = 3;
        % crop box
        cropBox = R(b).BoundingBox;
        % sensetivity for masking
        sen = SEN(s);
        singleCot(cFileList{e},I,M,cropBox,sen,spurAmount,snipAmount,skeletonDilateAmount,totalMaskErodeAmount)

        



            cropBox = R(b).BoundingBox;
            singleCot(I,M,cropBox);



        subI = imcrop(I,);
            subI = imcrop(I,R(b).BoundingBox);
            oI = subI;
            subI = imfilter(subI,fspecial('gaussian',[31 31],3));

            para.scales.value = 1;
            para.resize.value = 1;
            tmp = surKur(subI,para);
            subI = tmp(:,:,2);
            subI = imcomplement(subI);
            subI = bindVec(subI);

            %subI = imfilter(subI,fspecial('gaussian',[31 31],3));




            subM = imcrop(M,R(b).BoundingBox);
            subM = imfill(subM,'holes');




            cidx = find(subM==1);
            sample = subI(cidx);
            th = graythresh(sample);
            th = th + mag*th;
            sidx = find(sample > th);
            tmpM = zeros(size(subM));
            tmpM(cidx(sidx)) = 1;

            % .5 too sen
            sen = SEN(s);
            T = adaptthresh(subI,sen,'ForegroundPolarity','bright');
            tmpM = imbinarize(subI,T);


            tmpM = bwlarge(tmpM,2);
            tmpM = ~tmpM;
            tmpM = ~bwareaopen(tmpM,100);

            %tmpM = imfill(tmpM,'holes');




            skeleton = bwmorph(tmpM,'skeleton',inf);
            skeleton = bwmorph(skeleton,'spur',5);
            skeleton = bwmorph(tmpM,'skeleton',inf);


            bp = bwmorph(skeleton,'branchpoints');
            ep = bwmorph(skeleton,'endpoints');


            bPoints = [];
            ePoints = [];
            sPoints = [];
            path = {};
            pathcost = [];

            [bPoints(:,1),bPoints(:,2)] = find(bp);
            [ePoints(:,1),ePoints(:,2)] = find(ep);
            [sPoints(:,1),sPoints(:,2)] = find(skeleton);

            for v = 1:size(ePoints,1)
                stri = v;
                parfor w = 1:size(bPoints,1)
                    stpi = w;

                    q1 = find(sPoints(:,1) == ePoints(stri,1) & sPoints(:,2) == ePoints(stri,2));
                    q2 = find(sPoints(:,1) == bPoints(stpi,1) & sPoints(:,2) == bPoints(stpi,2));


                    Adj = Radjacency(sPoints',2^.5);
                    [path{v,w} , pathcost(v,w)] = dijkstra(Adj, q1 , q2);
                end
                v
            end
            [pcost,pidx] = min(pathcost,[],2);
            ridx = find(pcost < 7);

            for r = 1:numel(ridx)
                tmpX = sPoints(path{ridx(r),pidx(ridx(r))},2);
                tmpY = sPoints(path{ridx(r),pidx(ridx(r))},1);
                PTH = [tmpY tmpX];
                PTH = setdiff(PTH,bPoints,'rows');
                for p = 1:size(PTH,1)
                    skeleton(PTH(p,1),PTH(p,2)) = 0;
                end
                %plot(PTH(:,2),PTH(:,1),'k.')
            end






            skeleton = bwmorph(skeleton,'skeleton',inf);
            bp = bwmorph(skeleton,'branchpoints');
            ep = bwmorph(skeleton,'endpoints');

            bPoints = [];
            ePoints = [];
            sPoints = [];

            [bPoints(:,1),bPoints(:,2)] = find(bp);
            [ePoints(:,1),ePoints(:,2)] = find(ep);
            [sPoints(:,1),sPoints(:,2)] = find(skeleton);



            %{
            cLoop = ~skeleton;
            cLoop = cLoop.*subM;
            cLoop = imerode(cLoop,strel('disk',2,0));
            %}


            cLoop = ~imdilate(skeleton,strel('disk',3));
            useM = bwlarge(subM);
            for bor = 1:4
                useM(:,1) = 0;
                useM = imrotate(useM,90);
            end
            useM = imerode(useM,strel('disk',3,0));
            cLoop = cLoop.*useM;
            cR = regionprops(cLoop);
            cLoopLabel = bwlabel(cLoop);
            cLoopLabel = cLoopLabel + 1;
            resetIdx = find(useM == 0);
            cLoopLabel(resetIdx) = 0;
            rgbC = label2rgb(cLoopLabel);
            rgbC = double(rgbC)/255;





            out = flattenMaskOverlay(oI,logical(subM),.3,'r');
            out = flattenMaskOverlay(out,logical(tmpM),.3,'b');
            out = flattenMaskOverlay(out,logical(skeleton),.3,'y');
            oRGB = cat(3,oI,oI,oI);
            close all
            imshow([out rgbC oRGB],[]);
            hold on
            plot(bPoints(:,2),bPoints(:,1),'r.')
            plot(ePoints(:,2),ePoints(:,1),'b.')

            %{
            for r = 1:numel(ridx)
                tmpX = sPoints(path{ridx(r),pidx(ridx(r))},2);
                tmpY = sPoints(path{ridx(r),pidx(ridx(r))},1);
                PTH = [tmpY tmpX];
                PTH = setdiff(PTH,bPoints,'rows')
                plot(PTH(:,2),PTH(:,1),'k.')
            end
            %}

            drawnow


            snipFile = [pth filesep nm '_' num2str(b) '_' num2str(sen) '.jpg'];

            saveas(gca,snipFile);

            hold off
            %waitforbuttonpress


        end
        %out = flattenMaskOverlay(I,M);
        %imshow(out,[]);
        %drawnow
    end
end
