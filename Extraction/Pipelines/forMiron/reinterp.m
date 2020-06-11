function [V] = reinterp(F,X,Z)
    V = reshape(ba_interp2(F,X.E(1,:),X.E(2,:)),Z);
    funcF(V
end