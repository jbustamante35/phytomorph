function [m] = bwareaopenSmall(bw,areaSmall,conn)
    if nargin == 2;conn = 8;end
    CC = bwconncomp(bw,conn);
    R = regionprops(CC,'PixelIdxList','Area');
    kidx = ([R.Area] > areaSmall);
    R = R(kidx);
    m = zeros(size(bw));
    for e = 1:numel(R)
        m(R(e).PixelIdxList) = 1;
    end
end