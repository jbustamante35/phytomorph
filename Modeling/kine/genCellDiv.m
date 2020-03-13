function [cd] = genCellDiv(X,f,kw)
    cd = zeros(size(X));
    cd(1:f) = 1;
    cd = imfilter(cd,fspecial('gaussian',[1 kw(1)],kw(2)));
    
end