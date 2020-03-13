function [d] = makeRotatedCropBox(W,L,R,DX)
    [d1,d2] = ndgrid(linspace(-L/2,L/2,2),linspace(-W/2,W/2,2));
    d = [d1(:) d2(:) ones(size(d1(:)))];
    R = R + pi/2;
    T = [[cos(-R) sin(-R) DX(2)];[-sin(-R) cos(-R) DX(1)];[0 0 1]];
    d = T*d';
end