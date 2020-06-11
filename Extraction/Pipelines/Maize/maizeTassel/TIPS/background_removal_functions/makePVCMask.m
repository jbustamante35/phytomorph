function [pvc_mask] = makePVCMask(img)
%MAKEPVCMASK Summary of this function goes here
%   Detailed explanation goes here

si = std(img, 1, 3);
% Set super low variance pixels to the value of the 0.01 quantile
q = quantile(si(:), 0.01);
si_zero = si(:) <= q;
si(si_zero) = q;
mu = mean(img, 3);

tassel_norm = mu .* (si.^-1);
tassel_norm = bindVec(tassel_norm);

pvc_mask = tassel_norm > graythresh(tassel_norm);
pvc_mask = imopen(pvc_mask, strel('disk', 11, 0));

end

