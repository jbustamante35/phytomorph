function [f] = sampleApply(I,points,domain,func)
    if isempty(func)
        NZ_O = [size(domain,1),size(domain,2),size(domain,4)];
    else
        d = cat(3,domain(:,:,1,:) + points(1,1),domain(:,:,2,:) + points(1,2));
        tmp = squeeze(ba_interp2(I,d(:,:,1,:),d(:,:,2,:)));
        tmp = func(tmp);
        NZ_O = size(tmp);
    end
    f = zeros(prod(NZ_O),size(points,1));
    for pt = 1:size(points,1)
        d = cat(3,domain(:,:,1,:) + points(pt,1),domain(:,:,2,:) + points(pt,2));
        sample = squeeze(ba_interp2(I,d(:,:,1,:),d(:,:,2,:)));
        if ~isempty(func)
            sample = func(sample);
        end
        f(:,pt) = squeeze(sample(:));
    end
    f = reshape(f,[NZ_O size(points,1)]);
end






















