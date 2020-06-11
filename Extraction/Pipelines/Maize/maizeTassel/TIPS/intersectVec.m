function [] = intersectVec(v1,v2)
    options = optimset('MaxFunEvals',10000);
    f1 = @(a)(v1*[a*ones(2,1);1]);
    f2 = @(b)(v2*[b*ones(2,1);1]);
    func = @(x)sum(abs(f1(x(1)) - f2(x(2))));
    %x = fminunc(func,[1;1]);
    x = fminsearch(func,[1;1],options);
    p1 = f1(x(1));
    p2 = f2(x(2));
    P = [p1';p2'];
    P = mean(P,1);
end