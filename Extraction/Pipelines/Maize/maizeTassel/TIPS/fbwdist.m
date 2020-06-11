function [dt] = fbwdist(mask)
    here = 1;
    
    R = regionprops(mask);
    mask = imcrop(mask,R(1).BoundingBox);
    mask = imresize(mask,2);
    
    
    
    db = bwboundaries(mask,8,'holes');
    for e = 1:numel(db)
        [F,K] = smoothTasselContour(db{1},5);
        
    end
    
    
    z = zeros(size(mask));
    for e = 1:size(db{1},1)
        z(db{1}(e,1),db{1}(e,2)) = 1;
    end
    dt = bwdist(z);
    dt = double(mask).*dt - double(~mask).*dt;
    
    [dx,dy] = gradient(dt);
    
    
    
    
    
    %dt = imfilter(dt,fspecial('gaussian',[21 21],5),'replicate');
    
    mag = -1;
    for e = 1:20
        dt = dt + mag*del2(dt);
        %dt = imfilter(dt,fspecial('gaussian',[21 21],5),'replicate');
    end
    
    RGB = cat(3,dt>0,mask,zeros(size(mask)));
    
    iTH = 0;
    NO = inf;
    
    while NO ~= 1
        reg = dt > iTH;
        R = regionprops(reg);
        NO = numel(R);
        iTH = iTH + 1;
        %{
        imshow(reg,[]);
        drawnow
        waitforbuttonpress
        %}
    end
    
    outMask = dt > (iTH - 1);
end