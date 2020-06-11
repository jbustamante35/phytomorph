function [foreground] = schistosomaDetection(plateMask,background,I,minObjectSize)
    idx = find(plateMask);
    foreground = plateMask.*(I - background);
    v = foreground(idx) > graythresh(foreground(idx));
    foreground(idx) = v;
    foreground = bwareaopen(foreground,minObjectSize);
end