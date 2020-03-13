function [] = mainKineSim()
    close all
    
    % get a velocity function - default
    v = getVelocityFunction();
    
    
    % min/max length along the midline - from a position P to the apex
    minLength = 0;
    maxLength = 10;
    npLength = 301;
    % min/max width along the midline - from -1 to 1 as of 10.22.2019
    minWidth = -1;
    maxWidth =1;
    npWidth = 101;
    % render the grid
    rG = rGrid(minLength,maxLength,npLength,minWidth,maxWidth,npWidth,v);
   

    
    % render the velocity for a quick plot
    
    vv = rG.getVelocity();
    initState = rG.state;
    rG.stateT({'unfold','apex'});
   
    x = rG.grid(:,:,1);
    
    plot(x((end-1)/2,:),vv((end-1)/2,:),'k')
    hold on
    plot(x(1,:),vv(1,:),'r')
    plot(x(end,:),vv(end,:),'b')
    waitforbuttonpress
    close all
    h = figure;
    
    rG.stateT(initState);
    %rG.render(h);
    
    
    
    rG.updateGrid(30,.1,h,0);
    
    %{
    for loop = 1:30
        
        renderFluidGrid(X,sz,h,0,0);
        
        X = updatePosition(X,v,.1,sz);
        
        [l,w,midlineLength] = kineMetrics(X,sz);
    end
%}


end