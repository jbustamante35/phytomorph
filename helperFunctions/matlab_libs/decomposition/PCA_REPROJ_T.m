function [C ERR mu sd] = PCA_REPROJ_T(M,E,U,disp)
    if nargin == 3;disp = true;end
    ERR = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % reproject 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUTS:  
    %           M       : = data matrix
    %           E       : = basis vectors 
    %           U       : = mean of data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUTS:  
    %           C       : = unique finger print
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if disp;fprintf(['PCA_RE:start:taking off mean \n']);end
    M = bsxfun(@minus,M,U);
    if disp;fprintf(['PCA_RE:end:taking off mean \n']);end

    % subtract the mean
    %for i = 1:size(M,1)
    %    M(i,:) = M(i,:) - U;    % subtract the mean
    %end
    if disp;fprintf(['PCA_RE:start:proj into \n']);end
    % project to the "smaller" - (rotate) - subspace
    try
        C = mtimesx(E,'T',M);                    % project to get the coeffs
    catch
        C = E'*M;
    end
    if disp;fprintf(['PCA_RE:end:proj into \n']);end
    if nargout >= 2
        Mi = PCA_BKPROJ_T(C,E,U);
        ERR = sum((bsxfun(@plus,M,U) - Mi).^2,1).^.5;
    end
    if nargout > 2
        delta = (M - Mi);
        mu = mean(delta,1);
        sd = std(delta,1,1);
    end
end

