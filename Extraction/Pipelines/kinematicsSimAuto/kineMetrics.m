function [length,width,midlineLength] = kineMetrics(X,sz)

    X = reshape(X,sz);
    
    x = X(:,:,end-1);
    y = X(:,:,end);
    
    [dxdw,dxdl] = gradient(x);
    [dydw,dydl] = gradient(y);
    
    dw = (dxdw.^2 + dydw.^2).^.5;
    width = sum(dw,2);
    
    dl = (dxdl.^2 + dydl.^2).^.5;
    length = sum(dl,2);
    
    
    
    midlineLength = length((end-1)/2);
    
    
    
end