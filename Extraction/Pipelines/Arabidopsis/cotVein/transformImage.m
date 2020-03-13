function [ti] = transformImage(moving,fixed,x,trans,x0,nsz)
    
    % x0 - the reference point for the transformation
    % x - the transformation parameters

    if nargin == 5;nsz = size(moving);end
    if size(trans,2) == 5
        T = buildTrans(trans+x0);
    else
        T = trans;
    end

    x = mtimesx(T,x,'t');
    szT = size(fixed)/2;

    % move the points
    x(2,:) = x(2,:) + szT(2);
    x(1,:) = x(1,:) + szT(1);

    ti = ba_interp2(fixed,x(2,:),x(1,:));
    ti = reshape(ti,nsz);
end