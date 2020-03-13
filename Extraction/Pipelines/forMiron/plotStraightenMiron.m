function [] = plotStraighenMiron(X,funcS)
    func = funcS.func;
    dimX = funcS.dimX;
    dimY = funcS.dimY;

    X(dimY,:) = func(X(dimX,:));

end