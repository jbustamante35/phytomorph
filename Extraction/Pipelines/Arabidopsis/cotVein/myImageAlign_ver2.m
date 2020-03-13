function [x,df,fval] = myImageAlign_ver2(moving,movingMask,fixed,fixedMask,swarmSize,bounds,x,x0,disp,vwImage)

    




    [d] = generateAffineImageDomain(moving);
    [df] = generateAffineImageDomain(fixed);
    trans = [1 1 0 0 0];

    
    %edgeMoving = double(imdilate(edge(moving),strel('disk',3,0)));
    %edgeFixed = double(imdilate(edge(fixed),strel('disk',3,0)));
    %Rmoving = imresize(moving,.35);
    %Rmoving = moving;
    %edgeFixed = edge(fixed);
    %edgeMoving = edge(moving);
    %edgeMoving = edge(Rmoving);
   
    %[edgePointsM(:,1),edgePointsM(:,2)] = find(edgeMoving);
    %[edgePointsF(:,1),edgePointsF(:,2)] = find(edgeFixed);
    
    
    %edgePointsM = delaunayTriangulation(edgePointsM);

    %v = double(bwdist(moving < .8));



    func = @(x)myImageContrast_bwdistVersion(moving,fixed,d,x,x0,disp,edgePointsM,v);
    %func = @(x)myImageContrast(edgeMoving,edgeFixed,d,x);
    
    
    
    %func = @(x)myPointSetContrast(moving,fixed,d,x,edgePointsF);
    
    nvars = 5;
    
    options = optimoptions('particleswarm');
    options.Display = 'iter';
    options.UseParallel = true;
    options.UseVectorized = false;
    %options.HybridFcn = {@FMINCON};
    options.SwarmSize = swarmSize;
    
    x = particleswarm(func,nvars,b(1,:),b(2,:),options);
  
    options = optimset('fminsearch');
    options.Display = 'iter';
    x = fminsearch(func,x,options);

    options = optimoptions('fmincon');
    options.Display = 'iter';
    options.UseParallel = true;
    [x,fval] = fmincon(func,x,[],[],[],[],lb,ub,[],options);



    if disp
        imshow(fixed,[]);
        szf = size(fixed)/2;
        szm = size(moving)/2;
        hold on;
        dB = bwboundaries(moving > .8);
        T = buildTrans(x+x0);
        dB = dB{1};
        dB = bsxfun(@minus,dB,szm);
        dB = [dB ones(size(dB,1),1)];
        dB = (T*dB')';
        plot(x0(5)+szf(2),x0(4)+szf(1),'r*');
        plot(x0(5)+x(5)+szf(2),x0(4)+x(4)+szf(1),'g*');
        plot(dB(:,2)+szf(2),dB(:,1)+szf(1),'g')
        hold off
        drawnow
        here = 1;
    end
    %T = buildTrans(x);
    %T = buildTrans(x);
    %x = [1 1 0 0 0];
    %[fi] = transformImage(fixed,moving,df,inv(T));
    %imshowpair(fi, moving,'Scaling','joint');
    %imshowpair(fi, fixed,'Scaling','joint');
    
                    
end