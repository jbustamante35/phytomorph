function [selectMask] = generateSelectPoints(petriMask,pointSpacing)
    [p1,p2] = ndgrid(1:size(petriMask,1),1:size(petriMask,2));
    selectMask = (mod(p1,pointSpacing(1)) == 0) & (mod(p2,pointSpacing(2)) == 0);
    selectMask = petriMask .* selectMask;
end