function [p] = prob(I,density)
    szI = size(I);
    I = reshape(I,[prod(szI(1:end-1)), szI(end)]);
    x = density.binMapper(I);
    zidx = isnan(x);
    p = zeros(size(x));
    p(~zidx) = density.density(x(~zidx));
    p = reshape(p,szI(1:(end-1)));
end