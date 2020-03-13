function [data_out] = normalizeRawStack(data_in)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extract SIG from STACK
    % the signal that is processed is:
    % 1) subtract the initial signal from the stack
    % 2) divide by the max to have each signal eleOf [0,1]
    % 3) if the max == 0, then set to 1
    % 4) filter the signal
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % subtract each signal to start at zero
    data_out = bsxfun(@minus,data_in,data_in(:,:,1));
    % get the max value
    maxValue = max(data_out,[],3);
    % if max value == 0, set to 1
    maxValue(maxValue==0) = 1;
    % and have each signal go to one as the max
    data_out = bsxfun(@times,data_out,maxValue.^-1);
    sz = size(data_out);
    data_out = reshape(data_out,[prod(sz(1:2)) sz(3)]);
    % filter SIG
    data_out = imfilter(data_out,ones(1,11)/11,'replicate');
    % reshape SIG
    data_out = reshape(data_out,sz);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % renormalize
    % subtract filtered value
    data_out = bsxfun(@minus,data_out,data_out(:,:,1));
    % get the max value
    maxValue = max(data_out,[],3);
    % if max value == 0, set to 1
    maxValue(maxValue==0) = 1;
    % and have each signal go to one as the max
    data_out = bsxfun(@times,data_out,maxValue.^-1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end