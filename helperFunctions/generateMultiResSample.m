function [D] = generateMultiResSample(baseBoxSize,resolution,numberPoints)
    D = [];
    for r = 1:numel(resolution)
        baseBoxSizeTMP = baseBoxSize*resolution(r);
        dd1 = linspace(-baseBoxSizeTMP(1)/2,baseBoxSizeTMP(1)/2,numberPoints(1));
        dd2 = linspace(-baseBoxSizeTMP(2)/2,baseBoxSizeTMP(2)/2,numberPoints(2));
        [tD2 tD1] = ndgrid(dd2,dd1);
        D = [D;[tD1(:) tD2(:)]];
    end
end