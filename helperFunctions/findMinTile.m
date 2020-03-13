function [Z,M,indexImageMask,indexImage,LABEL] = findMinTile(IMR,stackMode,l1,l2,oldDX,toMatch,nl1,nl2,newDX,GR,functionBank,mfunctionBank,DEEP,Z,M,indexImageMask,indexImage)
    blendFunc = @(X)(1+exp(-.08*(X-5))).^-1;
    LABEL = zeros(size(indexImage));
    
    
    
    
    
    
    
    
    LABEL(find(indexImageMask)) = IMR(indexImage(find(indexImageMask)));
    for n = 1:size(newDX,1)
        oldMASK = zeros(size(Z));
        CP = zeros(size(oldMASK));
        
        nGR = bsxfun(@plus,GR,[nl1(newDX(n,1),newDX(n,2)) nl2(newDX(n,1),newDX(n,2))]);
       
        
        oldGR = {};
        for o = 1:size(oldDX,1)
            oldGR{o} = bsxfun(@plus,GR,[l1(oldDX(o,1),oldDX(o,2)) l2(oldDX(o,1),oldDX(o,2))]);
            IDXM = sub2ind(size(oldMASK),oldGR{o}(:,1),oldGR{o}(:,2));
            oldMASK(IDXM) = 1;
        end
        
        INNER_BLEND = bwdist(~oldMASK) - bwdist(oldMASK);
        INNER_BLEND = blendFunc(INNER_BLEND);
        OUTTER_BLEND = (1 - INNER_BLEND).*(oldMASK) + ~oldMASK;
        
        delta = 0;
        for e = 1:size(oldDX,1)
            tmpGR = oldGR{e};
            tmpMatch = toMatch(e);
            [iGR,ov1,ov2] = intersect(tmpGR,nGR,'rows');
            deltaN = bsxfun(@minus,functionBank(ov2,:),functionBank(ov1,tmpMatch));
            delta = delta + sum(deltaN.*deltaN,1).^.5;
        end
        
        [ssidx,sidx] = sort(delta);
        DEEPC = find(ssidx < 5);
        
        DEEP = min([max(DEEPC),DEEP]);
        DEEP
        
        if stackMode
            oldDX = [oldDX;newDX(n,:)];
            toMatch = [toMatch;sidx(DEEP)];
        end
            
        IDX = sub2ind(size(Z),nGR(:,1),nGR(:,2));
        newZ = NaN*zeros(size(Z));
        newM = NaN*zeros(size(M));
        newZ(IDX) = functionBank(:,toMatch(end));
        newM(IDX) = mfunctionBank(:,toMatch(end));
        Z = nanmean(cat(3,Z,newZ),3);
        M = nanmean(cat(3,M,newM),3);
        
        
        out = flattenMaskOverlay(Z,M > .8);
        imshow(out,[]);
        
        
        indexImageMask(newDX(n,1),newDX(n,2)) = 1;
        indexImage(newDX(n,1),newDX(n,2)) = sidx(DEEP);
        
        
        drawnow
        %waitforbuttonpress
    end

end