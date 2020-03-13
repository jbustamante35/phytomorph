function [X] = quantumEncode(MAX_B,INT,toSparse)
    % 0 - 255
    % the tensor product of the expanded form
    % of a binary vector in a truth table
    %--------------------------------------
    % b0 b1 | b0-ex b1-ex | dirc-form | N |
    %-------|-------------|-----------|---|
    % 0  0  | [0 1] [0 1] | [0 0 0 1] | 1 |
    % 0  1  | [0 1] [1 0] | [0 0 1 0] | 2 |
    % 1  0  | [1 0] [0 1] | [0 1 0 0] | 3 |
    % 1  1  | [1 1] [1 1] | [1 0 0 0] | 4 |
    %--------------------------------------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin == 2
        toSparse = false;
    end
    if ~toSparse
        BYTE_SIZE = 3;                              % log size of a byte
        MAX_WID = max(1,2^MAX_B / 2^BYTE_SIZE);     % number of bytes to allocate wide
        TALL = size(INT,1);
        X = zeros(TALL,MAX_WID,'uint8');
        % byte offset for index value into dirac-quantum form
        BYTE = floor(double(INT)*(2^BYTE_SIZE)^-1)+1;
        %BYTE = rem(double(INT),(2^MAX_B)+1)+1;
        % bit offset for index value into dirac-quantum form
        REM = rem(INT,2^BYTE_SIZE)+1;
        % position index value
        POS = sub2ind(size(X),(1:size(X,1))',BYTE);
        % set the bit
        X(POS) = bitset(X(POS),REM);
    else
        X = sparse(1:numel(INT),double(INT)+1,1);
    end
end