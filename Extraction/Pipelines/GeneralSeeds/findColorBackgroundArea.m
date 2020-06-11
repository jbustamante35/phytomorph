function [colorArea,cidx] = findColorBackgroundArea(I,darkT)
    LAB = rgb2lab(I);
    dark = LAB(:,:,1) < darkT;
    darkEdge = (imclearborder(dark) == 0) & (dark == 1);
    colorArea = ~darkEdge;
    cidx = find(colorArea);
end
