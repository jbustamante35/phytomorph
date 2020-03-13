function [F,X] = zoomOnPoint(func,P,W)
    try

        for dim = 1:size(P,2)
            xP{dim} = linspace(P(dim)-W.width(dim),P(dim)+W.width(dim),W.numP(dim));
        end
        X = cell(1,size(P,2));
        [X{:}] = ndgrid(xP{:});
        X = cell2mat(cellfun(@(X)X(:),X,'UniformOutput',false));
        F = func(X);
        F = reshape(F,[W.numP size(F,2)]);
    catch ME
        ME
    end
    
    
end