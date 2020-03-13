function [alignedPairMasks,alignedPairImage] = alignMasks_ver2(I,M,R)

        

        cropDims = 1.5*[[R.MajorAxisLength];[R.MinorAxisLength]]';
        cropLoc = reshape([R.Centroid],[2 numel(R)])';
        cropTheta = (pi/180)*[R.Orientation];
        


        aM = {};
        for b = 1:numel(R)
            aM{b} = imcrop(M,R(b).BoundingBox);
            
            
            
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
            imshow(aM{b},[]);
            drawnow
        end
        %{
        [optimizer, metric] = imregconfig('monomodal');
        optimizer.MaximumIterations = 100;
        
        [optimizer, metric] = imregconfig('multimodal');
        
        optimizer.InitialRadius = optimizer.InitialRadius*.01;
        optimizer.Epsilon = 1.5e-4;
        optimizer.GrowthFactor = 1 + (optimizer.GrowthFactor-1)*.5;
        optimizer.MaximumIterations = 100;
        %}
        
        dispO = false;
        
        sx = linspace(.8,1.2,5);
        sy = linspace(.8,1.2,5);
        
        
        
        
        PAD = 100;
        for template = 1:numel(aM)
            %{
            oTemplate = imcrop(I,R(template).BoundingBox);
            oTemplate = padarray(oTemplate,[PAD PAD],0,'both');
            %}
            
            mTemplate = padarray(aM{template},[PAD PAD],0,'both');
            oTemplate = padarray(aI{template},[PAD PAD],0,'both');
            RD = imref2d(size(mTemplate));
            
            
            for data = 1:numel(aM)
                MM = [];
                maxV = Inf;
                if template~= data
                    %oData = imcrop(I,R(data).BoundingBox);
                    %oData = padarray(oData,[PAD PAD],0,'both');
                    
                    mData = padarray(aM{data},[PAD PAD],0,'both');
                    oData = padarray(aI{data},[PAD PAD],0,'both');

                    %{
                    Trans_init = imregcorr(mData,mTemplate,'similarity');
                    mData = imwarp(mData,Trans_init,'OutputView',RD);
                    oData = imwarp(oData,Trans_init,'OutputView',RD);
                    imshowpair(oTemplate, overLay,'Scaling','joint');
                    %}
                    
                    %MoverLay = imregister(mData,mTemplate,'similarity',optimizer,metric,'DisplayOptimization',true,'PyramidLevels',4);
                    
                    [Trans,df] = myImageAlign(mData,mTemplate,100);
                    %{
                    Trans = imregtform(mData,mTemplate,'similarity',optimizer,metric,'DisplayOptimization',dispO,'PyramidLevels',4);
                    L = 1;
                    [~,~,Trans] = iterAlign(mData,mTemplate,oData,RD,L);
                    Trans = Trans{1};
                    overLay = imwarp(oData,Trans,'OutputView',RD);
                    MoverLay = imwarp(mData,Trans,'OutputView',RD);
                    %}
                    T = buildTrans(Trans);
                    [overLay] = transformImage(mTemplate,oData,df,inv(T));
                    [MoverLay] = transformImage(mTemplate,mData,df,inv(T));

                    %imshowpair(aM{template}, MoverLay,'Scaling','joint');
                    imshowpair(oTemplate, overLay,'Scaling','joint');
                    imshowpair(mTemplate, MoverLay,'Scaling','joint');
                    %waitforbuttonpress
                    alignedPairMasks{template,data} = cat(3,MoverLay,mTemplate);
                    alignedPairImage{template,data} = cat(3,overLay,oTemplate);
                    
                    %{
                    maxV = norm(MoverLay(:) - mTemplate(:));
                    alignedPair{template,data} = cat(3,MoverLay,mTemplate);
                    maxInit = maxV;
                    for x1 = 1:numel(sx)
                        for x2 = 1:numel(sy)
                            %{
                            dx = size(overLay)/2;
                            dx = flipdim(dx,2);
                            dx = dx - [sx(x1) sy(x2)].*dx;
                            tmpT = [[sx(x1) 0 0];[0 sy(x2) 0];[0 0 1]];
                            tmpT(3,1:2) = dx;
                            T = affine2d(tmpT);

                            %overLay2 = imwarp(overLay,T,'OutputView',RD);
                            %MoverLay2 = imwarp(MoverLay,T,'OutputView',RD);

                            overLay2 = imwarp(overLay,T,'OutputView',RD);
                            MoverLay2 = imwarp(MoverLay,T,'OutputView',RD);
                            
                            %}
                            
                            [overLay2,MoverLay2] = stretchImage(overLay,MoverLay,sx(x1),sy(x2),200,RD);
                            MoverLay2 = double(MoverLay2);
                            

                            %{
                            [optimizer, metric] = imregconfig('multimodal');
                            optimizer.InitialRadius = optimizer.InitialRadius*.1;
                            optimizer.Epsilon = 1.5e-4;
                            optimizer.GrowthFactor = 1 + (optimizer.GrowthFactor-1)*.5;
                            optimizer.MaximumIterations = 50;
                            
                            Trans2 = imregtform(MoverLay2,mTemplate,'similarity',optimizer,metric,'DisplayOptimization',dispO,'PyramidLevels',4);
                   
                            MoverLayFinal = imwarp(MoverLay2,Trans2,'OutputView',RD);
                            overLayFinal = imwarp(overLay2,Trans2,'OutputView',RD);
                            %}

                            [Trans2,df] = myImageAlign(MoverLay2,mTemplate);
                            T = buildTrans(Trans2);
                            [MoverLayFinal] = transformImage(mTemplate,MoverLay2,df,inv(T));
                            [overLayFinal] = transformImage(mTemplate,overLay2,df,inv(T));



                            tStore{x1,x2} = Trans2;
    
                            

                            %MM(x1,x2) = MoverLayFinal(:)'*mTemplate(:);
                            
                            MM(x1,x2) = norm(MoverLayFinal(:) - mTemplate(:));

                            if MM(x1,x2) < maxV
                                maxV = MM(x1,x2);
                                sv = [sx(x1) sy(x2)];
                                imshowpair(oTemplate, overLayFinal,'Scaling','joint');
                                drawnow
                            end
                            
                            
                            %imshowpair(oTemplate, MoverLayFinal,'Scaling','joint');
                            %waitforbuttonpress
                        end
                    end

                    %imshowpair(oTemplate, overLay,'Scaling','joint');
                    %waitforbuttonpress
                    
                    
                    close all
                    imshowpair(oTemplate, overLayFinal,'Scaling','joint');
                    drawnow
                    %waitforbuttonpress 
                    
                    if maxV ~= maxInit
                        alignedPair{template,data} = cat(3,MoverLayFinal,mTemplate);
                    end
                    %{
                    runner = cat(3,runner,alignedPair{template,data});
                    imshow(mean(runner,3),[]);
                    drawnow
                    pause(1);
                    %}
                    %}
                end
               
            end
        end
end


%{

    myImageAlign(mData,mTemplate)

    szD = size(mData)/2;
    xp = linspace(-szD(1),szD(1),size(mData,1));
    yp = linspace(-szD(2),szD(2),size(mData,2));
    [d1,d2] = ndgrid(xp,yp);
    d = [d1(:) d2(:) ones(size(d1(:)))];
    trans = [1 1 0 0 0];



    [fi] = transformImage(mData,mTemplate,d,trans);
    

%}