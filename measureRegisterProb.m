function [prob] = measureRegisterProb(reg)
    prob = sum((reg.*conj(reg)).^.5,2);
end