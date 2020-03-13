function [d] = myImageContrast_bwdistVersion(trans,trans0,moving,movingMask,fixed,fixedMask,d,disp,dispi,dispImage)
    

    if (dispi)% | (mod(dispVar,10) == 0))
        imshow(dispImage,[]);
        szf = size(fixed)/2;
        szm = size(moving)/2;
        hold on;
        dB = bwboundaries(moving > .8);
        T = buildTrans(trans+trans0);
        dB = dB{1};
        dB = bsxfun(@minus,dB,szm);
        dB = [dB ones(size(dB,1),1)];
        dB = (T*dB')';
        plot(trans0(5)+szf(2),trans0(4)+szf(1),'r*');
        plot(trans0(5)+trans(5)+szf(2),trans0(4)+trans(4)+szf(1),'g*');
        plot(dB(:,2)+szf(2),dB(:,1)+szf(1),'b')
        hold off
        drawnow
    end

    cookieCutter = false;
    edgeCompare = false;
    countNumBlack = false;
    DTcompare = true;

    fixedImageInMovingFrame = transformImage(moving,fixed,d,trans,trans0);
    

    if edgeCompare
        E = edge(fixedImageInMovingFrame);
        [ePoints(:,1),ePoints(:,2)] = find(E);
        [~,dist2] = DT.nearestNeighbor(ePoints);
        [d2] = pdist2(ePoints,DT.Points,'euclidean','Smallest',1);
        %[d1] = pdist2(DT.Points,ePoints,'euclidean','Smallest',1);
        d1 = dist2;
        edgeDistance = mean(d1) + mean(d2);
        edgeDistance = edgeDistance + std(d1).^2 + std(d2).^2;
    else
        dist2 = [0 0];
    end

    if cookieCutter
        %%%%%%%%%%%%%%%%
        % note: cutting then running the dist t
        fixedImageInMovingFrameC = fixedImageInMovingFrame.*movingMask;
        fixedImageInMovingFrameC = double(bwdist(~(fixedImageInMovingFrameC>.8)));
    end

    if countNumBlack
        % black pixel in mask count
        if ~isempty(movingMask)
            curMask = ~(fixedImageInMovingFrame > .8);
            numBlack = sum(curMask(:).*movingMask(:));
        end
    end

    if DTcompare
        % need to send in the whole mask and recalcuate the distance transform
        fixedImageInMovingFrame = double(bwdist(~(fixedImageInMovingFrame>.8)));
    end
    
   
    if cookieCutter
        %%%%%%%%%%%%%%%%
        % average?
        fixedImageInMovingFrame = .5*(fixedImageInMovingFrameC + fixedImageInMovingFrame);
    end


    if DTcompare
        if ~isempty(movingMask)
            fixedImageInMovingFrame = fixedImageInMovingFrame.*movingMask;
        end

        if ~isempty(fixedMask)
            % bring the fixed mask into the moving frame and apply to the
            % moving image
            fixedMaskInMovingFrame = transformImage(moving,fixed,d,trans,trans0);
            moving = moving.*fixedMaskInMovingFrame;
        end
        DTimage = fixedImageInMovingFrame(:) - moving(:);
        DTdistance = norm(DTimage);
    end
    
    ratio = trans(1) * trans(2).^-1;

    if ratio < .9
        rDistance = 1000;
    else
        rDistance = 0;
    end
    d = DTdistance + rDistance;
    %d = DTdistance + edgeDistance + rDistance;
    %dist1 = ;

    %distT = normpdf(dist1,0,3);
    %d = -sum(log(distT)) + numBlack + std(dist2).^2;
    %{
    fidx = find(moving ~= 0);
    [sidx] = dist1(fidx).^2;
    sidx = sort(sidx);
    toUse = round(perUse*numel(sidx));
    sidx = sidx(1:toUse);
    %}
    %d = norm(dist1);
    %d = norm(dist1);
    %d = norm(dist1) + mean(dist2) + std(dist2).^2 + numBlack;
    %d = norm(dist1);
    %d = -sum(log(distT));
   % d = sum(sidx).^.5; + mean(dist2);

    %itemplate = itemplate.*(vr > 0);

   % v = ba_interp2(itemplate,edgePointsM.Points(:,2),edgePointsM.Points(:,1));
    %d = sum(v(:));
    %d = norm(itemplate(:) - vr(:));
    %itemplate = imresize(itemplate,.35);
    %[templateP(:,1),templateP(:,2)] = find(edge(itemplate));
    %[~,delta1] = nearestNeighbor(edgePointsM,templateP);

    %tic
    %delta1 = pdist2(edgePointsM.Points,templateP,'euclidean','Smallest',1);
    %toc
    %d = mean(delta1);
    %delta2 = pdist2(templateP,edgePointsM.Points,'euclidean','Smallest',1);
    %d = mean(delta1) + mean(delta2) + extra;
    %d = bsxfun(@minus,itemplate,moving);
    %d = squeeze(sum(sum(d.*d,1),2).^.5);
    %d = norm(moving(:) - itemplate(:));
    %d = -moving(:)'*itemplate(:);
end