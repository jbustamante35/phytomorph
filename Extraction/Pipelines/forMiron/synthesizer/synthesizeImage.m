function [syn,count,BK,A] = synthesizeImage(backgroundImages,petriMasks,objects,maxOjects,N)
    squareSZ = [1000 1000];
    for e = 1:N
        e
        bki(e) = randi(size(backgroundImages,3),1);
        msk = petriMasks(:,:,bki(e));
        bk = backgroundImages(:,:,bki(e));
       
        

        R = regionprops(logical(msk),'BoundingBox');
        BK(:,:,e) = imresize(imcrop(bk,R(1).BoundingBox),squareSZ);
        tmpM = imresize(imcrop(msk,R(1).BoundingBox),squareSZ);
        
        midx = [];
        reducedM = bwdist(~msk) > 5;
        [midx(:,1),midx(:,2)] = find(reducedM==1);
        
        objectN(e) = randi(maxOjects,1);
        
        noiseN = 35;
        noiseM = .75;
        count(e) = 0;
        
        cp = zeros(size(bk));
        for o = 1:objectN(e)
            try
                
                oidx(o) = randi(numel(objects),1);
                objectI = objects(oidx(o)).I;
                pidx(o) = randi(size(midx,1),1);
                szI = size(objectI);
                strY = midx(pidx(o),1);
                stpY = strY + szI(1)-1;
                strX = midx(pidx(o),2);
                stpX = strX + szI(2)-1;
                sourceMask = double(objectI(:,:,3)>.5);
                targetMask = double(objectI(:,:,3)<.5);
                
                blendF = bwdist(sourceMask) < 2;
                blendIDX = find(~blendF);
                
                noiseIdx = randi(numel(blendIDX),noiseN,1);
                noiseValue = rand(noiseN,1);
                
                
                blendF(blendIDX(noiseIdx)) = noiseValue;
                
                sourceI = objectI(:,:,1);
                sourceCP = zeros(size(sourceI));
                sourceCP(round((size(sourceCP,1)-1)/2),round((size(sourceCP,2)-1)/2)) = 1;
                targetI = bk(strY:stpY,strX:stpX);
                
                resultI = targetI.*imcomplement(blendF) + sourceI.*blendF;
                
                cp(strY:stpY,strX:stpX) = sourceCP;
                bk(strY:stpY,strX:stpX) = resultI;
                %{
                imshow(bk,[]);hold on
                plot(strX,strY,'r.')
                waitforbuttonpress
                %}
                count(e) = count(e) + 1;
            catch ME
                ME
                
            end
        end
        cp = imdilate(cp,strel('disk',5,0));
        syn(:,:,e) = tmpM.*imresize(imcrop(bk,R(1).BoundingBox),squareSZ);
        A(:,:,e) = tmpM.*imresize(imcrop(cp,R(1).BoundingBox),squareSZ);
        syn(:,:,e) = tmpM*(syn(:,:,e) - BK(:,:,e));
        %imshow(bk,[]);
        %drawnow
        %waitforbuttonpress
    end
end