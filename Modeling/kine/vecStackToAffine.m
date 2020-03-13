function [a] = vecStackToAffine(tan,nor,delta)
    sz = size(tan');
    tan = reshape(tan',[sz(1) 1 sz(2)]);

    sz = size(nor');
    nor = reshape(nor',[sz(1) 1 sz(2)]);

    sz = size(delta');
    delta = reshape(delta',[sz(1) 1 sz(2)]);

    a = [tan,nor,delta];

    cap = zeros(1,size(a,2),size(a,3));
    cap(1,3,:) = 1;
    a = [a;cap];

    
end