function [tassel_mask] = makeTasselMaskKmeans(img)
%MAKETASSELMAKSKMEANS Summary of this function goes here
%   Detailed explanation goes here

szm = size(img);
tmp = reshape(img, [prod(szm(1:2)) szm(3)]);
kidx = kmeans(tmp, 2);
kidx = reshape(kidx, szm(1:2)) - 1; % subtract 1 to make binary

% Ensure that the tassel is ones
if sum(sum(kidx)) / numel(kidx) > 0.5
    kidx = -(kidx - 1);
end

tassel_mask = bwlarge(kidx);

end

