function [Y] = squeezeQbits2(X,bitSequence,LOG_BYTE_SIZE,TYPE)
    %{
    if numel(bitSequence) == 1
        bitSequence = bitSequence*ones(1,size(X,2));
    end
    %}
    bitNumber = 0;
    for e = 1:numel(bitSequence)
        bitNumber = bitNumber + numel(bitSequence{e});
    end


    NUMBER_BYTES = ceil(bitNumber/(2^LOG_BYTE_SIZE));
    Y = zeros(size(X,1),NUMBER_BYTES,TYPE);
    
    newColumn = 1;
    newBit = 1;
    for column = 1:size(X,2)
        for bit = 1:numel(bitSequence{column})
            bitValue = bitget(X(:,column),bitSequence{column}(bit));
            Y(:,newColumn) = bitset(Y(:,newColumn),newBit,bitValue);
            if newBit == (2^LOG_BYTE_SIZE)
                newColumn = newColumn + 1;
                newBit = 1;
            else
                newBit = newBit + 1;
            end
        end
    end
end