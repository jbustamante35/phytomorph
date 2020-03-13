function [F,theta] = getMidFrame(curve)
    midpoint = .5*(curve(end,:) - curve(1,:));
    T = midpoint / norm(midpoint);
    N = [T(2) -T(1)];
    midpoint = curve(1,:);
    F = [[T -T*midpoint'];[N -N*midpoint'];[0 0 1]];
    theta = atan2(T(1),T(2));
end