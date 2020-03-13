function [data_out] = rescaleData(data_in,N)
	sz = size(data_in);
    data_in = reshape(data_in,[prod(sz(1:2)) sz(3)]);
    X = (1:size(data_in,2))';
    Xi = linspace(1,size(data_in,2),N)';
    data_out = interp1(X,data_in',Xi);
    data_out = reshape(data_out',[sz(1:2) N]);
end