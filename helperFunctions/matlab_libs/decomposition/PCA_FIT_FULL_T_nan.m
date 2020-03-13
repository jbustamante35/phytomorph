function [S C U E L ERR LAM] = PCA_FIT_FULL_T_nan(M,COM)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % performs PCA (principal component analysis) using eigenvector
    % decomposition, backprojects to simulate the data and calculates the
    % error
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUTS:   
    %           M       : = data matrix
    %           COM     : = number of vectors 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUTS:  
    %           S       : = simulated signal (reconstruction)
    %           C       : = components ("unique" fingerprint)
    %           U       : = mean of data
    %           E       : = basis vectors 
    %           L       : = eig values
    %           ERR     : = error in reconstruction
    %           LAM     : = percent explained
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % take the mean
    fprintf(['PCA:start:taking off mean \n']);
    toOpDim = 2;
    U = nanmean(M,toOpDim);
    M = bsxfun(@minus,M,U);
    fprintf(['PCA:end:taking off mean \n']);
    % look at covariance
    fprintf(['PCA:start:creating COV \n']);
    COV = mynancov(M');
    fprintf(['PCA:end:creating COV \n']);
    % eig vector decomp
    fprintf(['PCA:start:decomposing COV \n'])
    [E L] = eigs(COV,COM);
    [J sidx] = sort(diag(L),'descend');
    E = E(:,sidx);
    L = diag(L);
    L = L(sidx);
    L = diag(L);
    fprintf(['PCA:end:decomposing COV \n'])
    % return eig values
    LAM = L;
    % calc percent explained
    L = cumsum(diag(L))*sum(diag(L))^-1;
    % get coeffs (fingerprints)
    for i = 1:size(E,2)
        for j = 1:size(M,2)
            C(i,j) = nansum(E(:,i).*M(:,j));
        end
    end
    % use the back-projection to create simulated signal
    S = PCA_BKPROJ_T(C,E,U);
    S(find(isnan(M))) = NaN;
    % calc the error
    ERR = nansum((S - M).^2,1).^.5;
end
