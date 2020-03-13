function [Y] = genCellD(X,para)
    Y = double(X < para(1));
    %Y = imfilter(Y,ones(1,para(2))/para(2),'replicate');
    %Y = imfilter(Y,ones(1,para(2))/para(2));
    Y = imfilter(Y,fspecial('gaussian',[1 2*para(2)],para(2)),'replicate');
    Y = para(3) * Y * max(Y(:))^-1;
end