function [n] = countHoles(holesImage)
    hR = regionprops(holesImage);
    n = numel(hR);
end