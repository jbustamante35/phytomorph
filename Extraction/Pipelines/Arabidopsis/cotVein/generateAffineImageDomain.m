function [d,sz] = generateAffineImageDomain(I,mag)
    if nargin == 1;mag = 1;end

    szD = size(I)/2;
    xp = linspace(-szD(1),szD(1),size(I,1)*mag);
    yp = linspace(-szD(2),szD(2),size(I,2)*mag);
    [d1,d2] = ndgrid(xp,yp);
    sz = size(d1);
    d = [d1(:) d2(:) ones(size(d1(:)))];

end