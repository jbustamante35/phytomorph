function [C,err] = PCA_REPROJ(M,E,U)
%%%%%%%%%%%%%%%%
% INPUTS:   M       : = data matrix
%           E       : = basis vectors 
%           U       : = mean of data
%%%%%%%%%%%%%%%%
% OUTPUTS:  C       : = unique finger print
%%%%%%%%%%%%%%%%
% subtract the mean
for i = 1:size(M,1)
    M(i,:) = M(i,:) - U;    % subtract the mean
end
% project to the "smaller" - (rotate) - subspace
%C = M*E;                    % project to get the coeffs
C = mtimesx(M,E);
if nargout == 2
    sim = PCA_BKPROJ(C,E,U);
    err = sim - M;
    err = sum(err.*err,2).^.5;
end
