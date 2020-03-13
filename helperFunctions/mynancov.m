function [cov] = mynancov(X)
    for i = 1:size(X,2)
        for j = 1:size(X,2)
            cov(i,j) = nanmean(X(:,i).*X(:,j));
        end
        i
    end
end