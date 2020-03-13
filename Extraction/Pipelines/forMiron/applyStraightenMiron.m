function [X] = applyStraightenMiron(X,funcS,combineFunction)
    for f = 1:numel(funcS)
        func = funcS(f).func;
        dimX = funcS(f).dimX;
        dimY = funcS(f).dimY;
        R(f,:) = func(X(dimX,:));
    end
    dataX = X(dimY,:);
    X(dimY,:) = combineFunction(dataX,R);
end