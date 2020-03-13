function [T] = transformFromGradient(globData,P)
    d(1) = ba_interp2(squeeze(globData(:,:,2,2)),P(1),P(2));
    d(2) = ba_interp2(squeeze(globData(:,:,3,2)),P(1),P(2));
    if norm(d) == 0
        d = rand(1,2);
    end
    d = d / norm(d);
    norD = [d(2) -d(1)];
    mkAffine = zeros(1,3);
    mkAffine(3) = 1;
    T = [[d',norD',P'];mkAffine];
end