function [X] = updatePosition_v2(X,V,dt,sz)
    v = V(X);
    
    
    
    dT = ones(size(v));
    dl = v;
    dw = zeros(size(dl));
    
    
    
    
    
    dX = v.*exp(-1j*X(:,4));
    
    
    
    
    
    
   
    
    vl = reshape(v,sz(1:2));
    vw = zeros(size(vl));
  
    
    
    
    [curlz,cav]= curl(vl,vw);
    
    
    
    
    
    %[dvdw,dvdl] = gradient(vl);
    %angVel = .5*dvdw;
    
    
    angVel = cav;
    da = angVel;
    %da = cumsum(angVel,1);
    da = da(:);
   
    dX = [dl dw dT da real(dX) imag(dX)]*dt;
   
    X = X + dX;
    
end