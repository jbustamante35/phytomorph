function [D] = generateMultiResSampleDomain_DISK(baseBoxSize,resolution,numberPoints)
    D = [];
    for r = 1:numel(resolution)
        baseBoxSizeTMP = baseBoxSize(1)*resolution(r);
        
        dd1 = linspace(0,baseBoxSizeTMP,numberPoints(1));
        dd2 = linspace(-pi,pi,numberPoints(2));
        [tD(:,:,2),tD(:,:,1)] = ndgrid(dd2,dd1);
        tD = cat(3,tD(:,:,1).*cos(tD(:,:,2)),tD(:,:,1).*sin(tD(:,:,2)));
        D = cat(4,D,tD);
    end
end