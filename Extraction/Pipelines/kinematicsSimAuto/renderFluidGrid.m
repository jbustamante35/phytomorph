function [] = renderFluidGrid(X,sz,h,theta,flag)



    X = reshape(X,sz);
    tmp = X(:,:,(end-1:end));
    tmpSZ = size(tmp);
    tmp = reshape(tmp,[prod(tmpSZ(1:2)) tmpSZ(3)]);
    tmp = [tmp ones(size(tmp,1),1)];
    
    p = squeeze(X(end,(end-1)/2,:));
    
    
    theta = X(end,(end-1)/2,4);
    
    offset = -pi/2;
    vT = [[cos(-theta) sin(-theta) p(end-1)];[sin(-theta) -cos(-theta) p(end)];[0 0 1]];
    aT = [[cos(theta+offset) sin(theta+offset) p(end-1)];[sin(theta+offset) -cos(theta+offset) p(end)];[0 0 1]];
    if flag
        aTa = inv(aT);
        tmp = mtimesx(aTa,tmp,'T')';
    end
   
    tmp = reshape(tmp(:,1:2),tmpSZ);
    
    
    
    Xgrid = reshape(X,sz);
    
    x = Xgrid(:,:,end-1);
    y = Xgrid(:,:,end);
    [dxdw,dxdl] = gradient(x);
    [dydw,dydl] = gradient(y);
    
    TAN_L = (dxdl.^2 + dydl.^2).^-.5;
    dxdl = dxdl.*TAN_L;
    dydl = dydl.*TAN_L;
    
    
    TAN_W = (dxdw.^2 + dydw.^2).^-.5;
    dxdw = dxdw.*TAN_W;
    dydw = dydw.*TAN_W;
    
    
    figure(h);
    for e = 1:size(tmp,1)
        plot(tmp(e,:,2),tmp(e,:,1),'c')
        hold on
    end
    plot(tmp(1,:,2),tmp(1,:,1),'b')
    
    
    for e = 1:size(tmp,2)
        plot(tmp(:,e,2),tmp(:,e,1),'m')
    end
    plot(tmp(:,1,2),tmp(:,1,1),'r')
    
    quiver(tmp(:,:,2),tmp(:,:,1),dydl,dxdl,'k')
    
    
    quiver(tmp(:,:,2),tmp(:,:,1),dydw,dxdw,'y')
    
    
    
    
    plot(vT(2,3),vT(1,3),'g.')
    quiver(vT(2,3),vT(1,3),vT(2,1),vT(1,1),10,'Color','g')
    
    %plot(0,0,'g*');
    
    %axis([-5 5 0 100]);
    hold off
    axis equal
    drawnow
    waitforbuttonpress
end