function [x,fval,df] = myImageAlign(x0,moving,movingMask,fixed,fixedMask,swarmSize,b,disp,dispi,dispBest,dispImage)
    
    
    boxSize = 350;
    xBK = x0;
    hfx = hsize(fixed);
    samplePoint = hfx + x0(4:5);
    cropBOX(1) = max(1,samplePoint(2)-boxSize);
    cropBOX(2) = max(1,samplePoint(1)-boxSize);
    cropBOX(3:4) = 2*boxSize;

    fixed = imcrop(fixed,cropBOX);
    dispImage = imcrop(dispImage,cropBOX);
    if ~isempty(fixedMask)
        fixedMask = imcrop(fixedMask,cropBOX);
    end

    % return the reference point to the image system
    x0(4:5) = x0(4:5) + hfx;
    x0(4:5) = x0(4:5) - flip((cropBOX(1:2) -1),2);
    x0(4:5) = x0(4:5) - hsize(fixed);
    

    
    %{
    imshow(fixed,[]);
    hold on;
    rectangle('Position',cropBOX);
    %}
   

    [d] = generateAffineImageDomain(moving);
    [df] = generateAffineImageDomain(fixed);
    %didx = find(movingMask > .8);

    %%%%%%%%%%%%%%%%%%
    % the edge method is not proving to be useful
    %{
    E = edge(movingMask);
    [ePoints(:,1),ePoints(:,2)] = find(E);
    DT = delaunayTriangulation(ePoints);
    %}
    DT = [];
    %%%%%%%%%%%%%%%%%%
    

    func = @(x)myImageContrast_ver2(x,x0,moving,movingMask,fixed,fixedMask,d,disp,dispi,dispImage);
    %func = @(x)myImageContrast(edgeMoving,edgeFixed,d,x);
    %func = @(x)myPointSetContrast(moving,fixed,d,x,edgePointsF);
    
    nvars = 5;
    dispFunc = @(X,Y)cotOptiOverlay(X,Y,[],dispImage,x0,fixed,moving,dispBest);
    options = optimoptions('particleswarm');
    if dispBest;options = optimoptions('particleswarm','OutputFcn',dispFunc);end
   
    options.Display = 'none';
    options.UseParallel = true & ~dispi;
    options.UseVectorized = false;
    options.FunctionTolerance = 1*10^-1;
    %options.HybridFcn = {@FMINCON};
    options.SwarmSize = swarmSize;
    lb = b(1,:);
    ub = b(2,:);
    [x,fval] = particleswarm(func,nvars,lb,ub,options);
    
    
    dispVar = 0;
    func = @(x)myImageContrast_bwdistVersion(x,x0,moving,movingMask,fixed,fixedMask,d,disp,dispi,dispImage);
    dispFunc = @(X,Y,Z)cotOptiOverlay(X,Y,Z,dispImage,x0,fixed,moving,dispBest,curveBank);
    options = optimset('OutputFcn',dispFunc);
    options = optimset('fminsearch');
    options.Display = 'none';
    %[x,fval] = fminsearch_phytoMorph(func,x,options);
    [x,fval] = fminsearch(func,x,options);
    

    %{
    options = optimoptions('fmincon');
    options.Display = 'iter';
    options.UseParallel =  true & ~dispi;
    [x,fval] = fmincon(func,x,[],[],[],[],lb,ub,[],options);
    %}
    %{
    rotVec = pi/4;
    dispFunc = @(X,Y,Z)cotOptiOverlay(X,Y,Z,dispImage,x0,fixed,moving,dispBest,curveBank);
    options = optimset('OutputFcn',dispFunc);
    options.Display = 'iter';
    xinit = x;
    for r = 1:4
        x = xinit;
        x(3) = x(3) + (r-1)*pi/2;
        xrot(r,:) = fminsearch_phytoMorph(func,x,options);
    end
    %}


    if disp
        imshow(dispImage,[]);
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
    end
    %T = buildTrans(x);
    %T = buildTrans(x);
    %x = [1 1 0 0 0];
    %[fi] = transformImage(fixed,moving,df,inv(T));
    %imshowpair(fi, moving,'Scaling','joint');
    %imshowpair(fi, fixed,'Scaling','joint');
    
    %x(4:5) = x(4:5) + hsize(fixed);
    %x(4:5) = x(4:5) + flip((cropBOX(1:2) -1),2);
                    
end