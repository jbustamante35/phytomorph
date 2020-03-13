function [S C U E L ERR LAM] = PCA_FIT_FULL(M,COM,verbose)
if nargin == 2;verbose=true;end
%%%%%%%%%%%%%%%%
% INPUTS:   M       : = data matrix
%           COM     : = number of vectors 
%%%%%%%%%%%%%%%%
% OUTPUTS:  S       : = simulated signal (reconstruction)
%           C       : = components ("unique" fingerprint)
%           U       : = mean of data
%           E       : = basis vectors 
%           L       : = eig values
%           ERR     : = error in reconstruction
%           LAM     : = percent explained
%%%%%%%%%%%%%%%%

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
    if verbose;fprintf(['PCA:start:taking off mean \n']);end
    toOpDim = 1;
    U = mean(M,toOpDim);
    M = bsxfun(@minus,M,U);
    if verbose;fprintf(['PCA:end:taking off mean \n']);end
    % look at covariance
    if verbose;fprintf(['PCA:start:creating COV \n']);end
    try
        COV = mtimesx(M,'T',M,'speed');
    catch
        COV = M'*M;
    end
    COV = COV / size(M,toOpDim);
    if verbose;fprintf(['PCA:end:creating COV \n']);end
    % eig vector decomp
    if verbose;fprintf(['PCA:start:decomposing COV \n']);end
    [E L] = eigs(COV,COM);
    
    
    [J sidx] = sort(diag(L),'descend');
    E = E(:,sidx);
    L = diag(L);
    L = L(sidx);
    L = diag(L);
    
    if verbose;fprintf(['PCA:end:decomposing COV \n']);end
    % return eig values
    LAM = L;
    % calc percent explained
    L = cumsum(diag(L))*sum(diag(L))^-1;
    % get coeffs (fingerprints)
    try
        C = mtimesx(M,E);
    catch
        C = M*E;
    end
    % use the back-projection to create simulated signal
    S = PCA_BKPROJ(C,E,U);
    % calc the error
    ERR = sum((S - M).^2,1).^.5;


%{
% take the mean
fprintf(['PCA:start:taking off mean \n']);
U = mean(M,1);
for i = 1:size(M,1)
    M(i,:) = M(i,:) - U;
end
fprintf(['PCA:end:taking off mean \n']);
% look at covariance
COV = cov(M);
% eig vector decomp
[E L] = eigs(COV,COM);
% return eig values
LAM = L;
% calc percent explained
L = cumsum(diag(L))*sum(diag(L))^-1;
% get coeffs (fingerprints)
if isdeployed
    C = M*E;
else
    C = mtimesx(M,E);
end
% use the back-projection to create simulated signal
S = PCA_BKPROJ(C,E,U);
% calc the error
ERR = sum((S - M).^2,2).^.5;
    %}