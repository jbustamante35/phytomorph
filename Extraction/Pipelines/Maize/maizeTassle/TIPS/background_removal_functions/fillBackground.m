function [img] = fillBackground(srcMask,targetMask, img)
%FILLBACKGROUND Summary of this function goes here
%   Detailed explanation goes here

fill = img;
        for k = 1:size(img, 3)
            tmp = img(:,:,k);
            mu_k(k) = mean(tmp(find(~srcMask)));
            tmp = fill(:,:,k);
            tmp(find(targetMask)) = mu_k(k);
            fill(:,:,k) = tmp;
        end
        img = fill;
        
end

