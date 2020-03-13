function [m] = bwareaopenRange(bw,areaRange)
    R = regionprops(bw,'PixelIdxList','Area');
    kidx = ([R.Area] > areaRange(1) &  [R.Area] < areaRange(2));
    R = R(kidx);
    m = zeros(size(bw));
    for e = 1:numel(R)
        m(R(e).PixelIdxList) = 1;
    end
end