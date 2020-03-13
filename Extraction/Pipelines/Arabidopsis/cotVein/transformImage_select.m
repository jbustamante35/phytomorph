function [ti] = transformImage_select(moving,fixed,d,didx,x,x0,nsz)
    
    % x0 - the reference point for the transformation
    % x - the transformation parameters

    if nargin == 6;nsz = size(moving);end

    if size(x,2) == 5
        T = buildTrans(x+x0);
    else
        T = x;
    end

    d = mtimesx(T,d,'t');
    szT = size(fixed)/2;

    % move the points
    d(2,:) = d(2,:) + szT(2);
    d(1,:) = d(1,:) + szT(1);
    ti = zeros(size(fixed));

    ti = zeros(size(moving));

    ti(didx) = ba_interp2(fixed,d(2,:),d(1,:));
    
end