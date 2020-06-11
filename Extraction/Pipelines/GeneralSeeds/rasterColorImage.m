function [I] = rasterColorImage(I,func,cidx,blankV)
    szI = size(I);
    I = reshape(I,[prod(szI(1:2)) 1 szI(3)]);
    if nargin > 1
        data = squeeze(I(cidx,:,:));
        data = func(data);
        if size(data,2) ~= szI(end)
            udata = mean(data);
            I = blankV*ones(size(I,1),size(data,2));
            szI(end) = size(data,2);
        end
        I(cidx,:,:) = data;
        I = reshape(I,szI);
    end
end