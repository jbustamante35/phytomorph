function [PROGRAM_TEST] = generateInitProgram(X,MAX_P_NUMBER,N,MAX_WID)
    if issparse(X)
        PROGRAM_TEST = sparse(1,size(X,2));
        numBitToSet = randi(size(PROGRAM_TEST,2),1,1);
        bitsToSet = randi(size(PROGRAM_TEST,2),numBitToSet,1);
        for b = 1:numBitToSet
            PROGRAM_TEST(bitsToSet(b)) = 1;
        end
    else
        PROGRAM_TEST = randi(MAX_P_NUMBER,N,MAX_WID,'uint8');
    end
end