function [c,ceq] = abc(noncon,x)
    c = [];
    ceq = noncon(x);
end