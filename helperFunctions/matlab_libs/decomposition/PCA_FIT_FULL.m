function [S C U E L ERR LAM] = PCA_FIT_FULL(M,COM)
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
    U = mean(M,1);
    M = bsxfun(@minus,M,U);
    %{
    for i = 1:size(M,1)
        M(i,:) = M(i,:) - U;
    end
    %}
    % look at covariance
    COV = cov(M);
    %[COV,sigma] = corrcov(COV);
    % eig vector decomp
    [E L] = eigs(COV,COM);
    % return eig values
    LAM = L;
    % calc percent explained
    L = cumsum(diag(L))*sum(diag(L))^-1;
    % get coeffs (fingerprints)
    C = M*E;
    % use the back-projection to create simulated signal
    S = PCA_BKPROJ(C,E,U);
    % calc the error
    ERR = sum((S - M).^2,2).^.5;
end