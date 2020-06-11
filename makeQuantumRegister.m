function [reg] = makeQuantumRegister(v)
    Q = numel(v);
    v = reshape(v,[Q/3 3]);
    n = size(v,1);
    c = [v(:,1).*exp(1i*v(:,2)) (1-v(:,1)).*exp(1i*v(:,3))]; 
    
    %{
    n = sum(c.*conj(c),2);
    c = bsxfun(@times,c,n.^-.5);
    n = sum(c.*conj(c),2);
    %}

    % init the register
    reg = c(1,:)';
    for e = 2:size(c,1)
        reg = reg*c(e,:);
        reg = reg(:);
    end
    
    %n = sum((reg.*conj(reg)).^.5,2);
    % measure eigenvalues
    %plot(n);
end