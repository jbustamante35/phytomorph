function [T] = buildTrans(trans)
    %z = zeros([1 1 size(trans,3)]);
    %S = [[trans(1,:,:) z z];[z trans(2,:,:) z];[z z ones(size(z)]];
    S = [[trans(1) 0 0];[0 trans(2) 0];[0 0 1]];
    th = trans(3);
    R = [[cos(th) sin(th) 0];[-sin(th) cos(th) 0];[0 0 1]];
    D = [[1 0 0];[0 1 0];[0 0 1]];
    D(1:2,3) = trans(4:5);
    T = D*R*S;
    %T = R;
end