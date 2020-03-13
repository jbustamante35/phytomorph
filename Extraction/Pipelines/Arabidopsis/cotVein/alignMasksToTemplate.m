function [dataPoint] = alignMasksToTemplate(fileName,I,M,R,mTemplate,N)

        

        cropDims = 1.5*[[R.MajorAxisLength];[R.MinorAxisLength]]';
        cropLoc = reshape([R.Centroid],[2 numel(R)])';
        cropTheta = (pi/180)*[R.Orientation];
        
        aM = {};
        for b = 1:numel(R)
            aM{b} = imcrop(M,R(b).BoundingBox);
            
            cropFunc{b} = @(I)imcrop(I,cropDims(b,2),cropDims(b,1),cropTheta(b),cropLoc(b,1:2));
            
            tmp = rimcrop(M,cropDims(b,2),cropDims(b,1),cropTheta(b),cropLoc(b,1:2));
            tmpI = rimcrop(I,cropDims(b,2),cropDims(b,1),cropTheta(b),cropLoc(b,1:2));
            
            tmp = tmp > .9;
            tmp = imclearborder(tmp);
            
            tmpI = tmpI.*tmp;
            
            
            E = edge(tmp);
            E = imdilate(E,strel('disk',5,0));
            hood = strel('disk',5,0);
            tmp = bwlarge(imfill(tmp > .9,'holes'));
            
            dB =  bwboundaries(tmp);
            dB = dB{1};
            para{1} = 15;
            out = cwtK_closed_imfilter(dB,para);
            sidx = sub2ind(size(tmp),dB(:,1),dB(:,2));
            ttmp = zeros(size(tmp));
            ttmp(sidx) = out.K;
            SM = sum(abs(ttmp),2);
            
            
            
            topH = SM(1:round((size(SM,1)/2)));
            botH = SM(round((size(SM,1)/2)):end);
            botH = botH / sum(botH);
            topH = topH / sum(topH);
            topH = flip(topH,1);
            
            TOP = topH'*(1:numel(topH))';
            BOT = botH'*(1:numel(botH))';
            
            %if sum(botH) > sum(topH)
            
            TOP - BOT;
            if TOP > BOT
                tmp = flip(tmp,1);
                tmpI = flip(tmpI,1);
            end
            
            aM{b} = double(tmp);
            aI{b} = double(tmpI);
            %imshow(aM{b},[]);
            %drawnow
        end
        
        
        


        PAD = 100;
        for data = 1:numel(aM)

            mData = padarray(aM{data},[PAD PAD],0,'both');
            oData = padarray(aI{data},[PAD PAD],0,'both');

            [Trans,df] = myImageAlign(mData,mTemplate,N);
            
            T = buildTrans(Trans);
            [overLay] = transformImage(mTemplate,oData,df,inv(T));
            [MoverLay] = transformImage(mTemplate,mData,df,inv(T));
            

            dataPoint(data).alignedMask = MoverLay;
            dataPoint(data).alignedImage = overLay;
            dataPoint(data).T = T;
            dataPoint(data).Mask = mData;
            dataPoint(data).Image = oData;
            dataPoint(data).cropFunc = cropFunc{data};
            dataPoint(data).fileName = fileName;

            imshowpair(mTemplate, MoverLay,'Scaling','joint');
            
        end
    here =1;
end
