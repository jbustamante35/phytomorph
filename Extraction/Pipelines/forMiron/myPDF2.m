function [Z] = myPDF2(data,N)
    Z = zeros(N);
    [Yx,binsX] = discretize(data(1,:),N(2));
    [Yy,binsY] = discretize(data(2,:),N(1));
    rdata = round(data);
    idx = sub2ind(size(Z),Yy,Yx);
    for e = 1:numel(idx)
        Z(idx(e)) = Z(idx(e)) + 1;
    end
    Z = Z / sum(Z(:));
end
