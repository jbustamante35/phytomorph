function [ret] = getVelocityFrame(I,x)
    v = mean(I,1);
    v = v(2:3);
    t = v / norm(v);
    n = [t(2) -t(1)];
    T = eye(3);
    T(:,1) = [t';0];
    T(:,2) = [n';0];
    
    
    
    ret = funcK.affine(x.P,x.globData,x.funcT,T);
end