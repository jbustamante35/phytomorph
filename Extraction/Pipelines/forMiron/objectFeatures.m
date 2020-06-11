function [objectList] = objectFeatures(I,vI,oM,vM,objectList,domain)
    div = divergence(vI(:,:,1),vI(:,:,2));
    cur = curl(vI(:,:,1),vI(:,:,2));
    oM = double(oM);
    u = vI(:,:,1);
    v = vI(:,:,2);



    domainCP = (domain.sz-1)/2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % work up on velocity object
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for e = 1:numel(objectList)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make a static reference frame, sample 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        TAN = [sin(objectList(e).Orientation*pi/180);cos(objectList(e).Orientation*pi/180)];
        NOR = [TAN(2);-TAN(1)];
        DX = [objectList(e).Centroid]';
        F = [TAN NOR DX];
        %F = [NOR TAN DX];
        samDomain = (F*domain.x')';
       	tmpSample = ba_interp2(oM,samDomain(:,1),samDomain(:,2));
        tmpSample = reshape(tmpSample,domain.sz);
        tmpSampleMask = tmpSample > .5;
        tmpSampleMask = imfill(~tmpSampleMask,domainCP+1,8) & tmpSampleMask;
        % balance the left and right
        leftH = tmpSampleMask(:,1:(end-1)/2);
        rightH = tmpSampleMask(:,((end-1)/2+2):end);
        lH = std(find(sum(leftH,2)));
        rH = std(find(sum(rightH,2)));
        % flip TAN if needed
        if lH > rH;F = [-TAN NOR DX];end
        % store frame and sample
        objectList(e).static_frame = F;
        samDomain = (objectList(e).static_frame*domain.x')';
        tmpSample = ba_interp2(oM,samDomain(:,1),samDomain(:,2));
        tmpSample = reshape(tmpSample,domain.sz);
        tmpSampleMask = tmpSample > .5;
        tmpSampleMask = imfill(~tmpSampleMask,domainCP+1,8) & tmpSampleMask;
        %tmpSample = tmpSample.* tmpSampleMask;

        % store samples
        objectList(e).static_samMask = tmpSampleMask(:);
        objectList(e).static_sam_I = tmpSample(:);
        objectList(e).static_sam_vx = ba_interp2(u,samDomain(:,1),samDomain(:,2));
        objectList(e).static_sam_vy = ba_interp2(v,samDomain(:,1),samDomain(:,2));
        %{
        imshow(I,[]);
        hold on
        imshow(tmpSample,[]);
        waitforbuttonpress
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get div and curl
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        objectList(e).div = mean(div(objectList(e).PixelIdxList));
        objectList(e).curl = mean(cur(objectList(e).PixelIdxList));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sample the centroid velocity
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        objectList(e).centroid_vx = ba_interp2(vI(:,:,1),objectList(e).Centroid(1),objectList(e).Centroid(2));
        objectList(e).centroid_vy = ba_interp2(vI(:,:,2),objectList(e).Centroid(1),objectList(e).Centroid(2));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make a velocity reference frame
        % and measure
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        TAN = [objectList(e).centroid_vx objectList(e).centroid_vy];
        TAN = TAN / norm(TAN);
        NOR = [TAN(2) -TAN(1)];
        % get the vel field in the blob
        tmpVelField = [u(objectList(e).PixelIdxList),v(objectList(e).PixelIdxList)];
        % get velocities relative to the frame
        alongDir = tmpVelField*TAN';
        alongNor = tmpVelField*NOR';
        % measure the expected velocities - (x,y)
        objectList(e).expected_vx = mean(tmpVelField(:,1));
        objectList(e).expected_vy = mean(tmpVelField(:,2));
        % measure the expected velocities - (u,v)
        objectList(e).expected_vt = mean(alongDir);
        objectList(e).expected_vn = mean(alongNor);
        % measure the std velocities - (u,v)
        objectList(e).var_vt = std(alongDir);
        objectList(e).var_vn = std(alongNor);
        % measure the centroid
        objectList(e).x = objectList(e).Centroid(1);
        objectList(e).y = objectList(e).Centroid(2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make a velocity reference frame
        % and measure
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        objectList(e).dynamic_frame = [[TAN' NOR' objectList(e).Centroid']];
        samDomain = (objectList(e).dynamic_frame*domain.x')';
        objectList(e).dynamic_sam_I = ba_interp2(I,samDomain(:,1),samDomain(:,2));
        objectList(e).dynamic_sam_vx = ba_interp2(u,samDomain(:,1),samDomain(:,2));
        objectList(e).dynamic_sam_vy = ba_interp2(v,samDomain(:,1),samDomain(:,2));
        objectList(e).domainSZ = domain.sz;
        % transpose data for use
        objectList(e).BoundingBox = objectList(e).BoundingBox';
        objectList(e).Centroid = objectList(e).Centroid';
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    here =1 ;
end
