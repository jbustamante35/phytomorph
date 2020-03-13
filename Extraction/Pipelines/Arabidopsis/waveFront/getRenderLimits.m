function [mM] = getRenderLimits(data)
    for l = 1:size(data,3)
        tmp = round(data(:,:,l));
        mM(l,1) = 1;
        mM(l,2) = max(tmp(:));
    end
end