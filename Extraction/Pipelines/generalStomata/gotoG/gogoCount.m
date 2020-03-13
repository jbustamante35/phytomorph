function [res,CNT] = gogoCount(nptSZ,fileName,npt,G,tE,tU,trainedNet,tree,GLM,STEP,beta,nb,pnet,FDA,fIDX,fIDX_step,disp,rPath)

    [pth,nm,ext] = fileparts(fileName);

    res = zeros(size(npt,1),1);
    I = imread(fileName);
   
    
    parfor p = 1:size(npt,1)
        R = zeros(1,1);
        if disp
            imshow(I,[]);
            hold on
            plot(npt(p,2),npt(p,1),'r*')
        end
        
        
        tmpD = cat(3,G(:,:,1,:)+npt(p,1),G(:,:,2,:)+npt(p,2));
        subI = ba_interp2(I,squeeze(tmpD(:,:,2,:)),squeeze(tmpD(:,:,1,:)));
        %subIN = bsxfun(@minus,subI,mean(subI,1));
        %subI = abs(fft(subI,[],1));
        %subI = subI(1:40,:,:);
        
        
        %{
        % complex work
        tmpD = cat(7,G(:,:,:,:,:,:,1)+npt(p,1),G(:,:,:,:,:,:,2)+npt(p,2));
        subI = ba_interp2(I,tmpD(:,:,:,:,:,:,2),tmpD(:,:,:,:,:,:,1));
        subI = abs(fft(subI,[],1));
        subI = subI(1:30,:,:,:);
        tmpSZ = size(subI);
        subI = reshape(subI,[tmpSZ(1:2) prod(tmpSZ(3:4))]);
        %}
        
        %subIN = abs(fft(subIN,[],1));
        %subIN = subIN(1:40,:,:);
        %[tC_tmp,tmpErr] = PCA_REPROJ(subIN(:)',tE,tU);
        
        %[tC_tmp] = PCA_REPROJ(subI(:)',tE,tU);
        %[tC_tmp,tmpErr] = PCA_REPROJ(subI(:)',tE,tU);
        %tC_tmp = [tC_tmp tmpErr];
        
        
        
        %class(GLM)
        
        %STACK_whole(:,:,:,p) = subI;
        %STACK_vec(p,:) = tC_tmp;
        
        
        
        % for CNN
        tmp = trainedNet.predict(subI);
        R(1) = tmp(2);
        
        
        %{
        R(2) = tree.predict(tC_tmp(fIDX));
        
        R(3) = GLM.predict(tC_tmp(fIDX));
        
        R(4) = STEP.predict(tC_tmp(fIDX_step));
        
        R(5) = [1 tC_tmp(:)']*beta;
        
        [~,tmp] = nb.predict(tC_tmp(fIDX));
        R(6) = tmp(2);
        
        
        
        
        %tmp = sim(pnet,tC_tmp(fIDX)');
        %R(7) = tmp(1);
        
        [~,tmp] = FDA.predict(tC_tmp(fIDX));
        R(8) = tmp(2);
        %}
        
        
        if disp
            for m = 1:size(subI,3)
                tmp = subI(:,:,m);
                tmp = bindVec(tmp);
                subI(:,:,m) = tmp;
            end

            imshow(subI,[]);
            hold off
            drawnow
        end
        
       
        
        %fprintf([num2str(p) ':' num2str(size(npt,1)) '\n']);
        res(p,:) = R;
    end
    
    subI = ba_interp2(I,npt(:,2),npt(:,1));
    subI = reshape(subI,nptSZ(1:2));
    sig = reshape(res,nptSZ(1:2));
    

    %sig = imfilter(sig,fspecial('gaussian',[71 71]),'replicate');
    %sig = imfilter(sig,fspecial('gaussian',[MM(m)]),'replicate');

    %sig = imfilter(sig,fspecial('gaussian',[3]),'replicate');
    sig = imfilter(sig,fspecial('gaussian',[5]),'replicate');
    %mm = sig > TH(l);
    %mm = sig > 0.53;
    mm = sig > .21111;
    %mm = sig > 0.1556;
    mm = bwareaopen(mm,40);
    
    
    
    
    
    
    R = regionprops(logical(mm),'Centroid','Area','Perimeter','PixelIdxList');


    mm = zeros(size(mm));
    mm2 = zeros(size(mm));
    fidx = count([R.Area])==1;
    fidx2 = count([R.Area])==2;

    
    
    
    
    
    
    fidx = find(fidx == 1);
    for f = 1:numel(fidx)
        mm(R(fidx(f)).PixelIdxList) = 1;
    end
    
    fidx2 = find(fidx2 == 1);
    for f = 1:numel(fidx2)
        mm2(R(fidx2(f)).PixelIdxList) = 1;
    end

    mm = logical(mm);
    mm2 = logical(mm2);
    
    
    CNT = numel(fidx) + 2*numel(fidx2);
    
    
    out = flattenMaskOverlay(bindVec(imresize(subI,size(mm))/255),mm,.55,'r');
    out = flattenMaskOverlay(out,mm2,.55,'g');
    
    
    mkdir('./output/');
    fileList{1} = ['./output/' nm '.jpg'];
    imwrite(out,fileList{1});
    fileList{2} = ['./output/' nm '_count.csv'];
    csvwrite(fileList{2},CNT);
    R = regionprops('table',logical(mm+mm2),'Centroid','Area','Perimeter','PixelIdxList');
    fileList{3} = ['./output/' nm '_data.csv']
    writetable(R,fileList{3})
    fileList{4} = ['./output/' nm '_raw.jpg'];
    imwrite(bindVec(imresize(subI,size(mm))/255),fileList{4});
    
    
    pushToiRods(rPath,fileList);
    
    
end