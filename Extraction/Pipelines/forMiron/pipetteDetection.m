function [pipetteArea] = pipetteDetection(plateMask,background,I,minObjectSize)
    idx = find(plateMask);
    foreground = plateMask.*(I - background);
    v = foreground(idx) > graythresh(foreground(idx));
    foreground(idx) = v;
    foreground = bwlarge(logical(foreground));
    foreground = bwareaopen(foreground,minObjectSize);
    pipetteArea = sum(foreground(:));
end