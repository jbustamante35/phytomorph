function [Y] = evalProgram(programBank,bitRegister)


    % inputs:- programBank -> a sequence of a set of programs to evaluate
    % br:- bitRegister -> a sequence of a set of bits [set,depth,sequence]

    % get the program sequence length
    programSequenceLength = size(programBank,3);
    registerSequenceLegnth = bitRegister.getSeqLength();

    programPointer = 1;
    registerPointer = 1;
    for registerSequence = 1:registerSequenceLegnth


        X = bitRegister.entangleRegisterAtN(registerPointer);


        % splice off Cout bits
        if registerSequence >= 2
            here = 1;
            X = bitRegister.entangleNewBits(registerPointer,Y(:,2:end,registerPointer-1));
        end



        for programNumber = 1:size(programBank,1)


            if ~issparse(X)
                Y(:,programNumber,registerPointer) = boo(X,programBank(programNumber,:,programPointer));
            else
                Y(:,programNumber,registerPointer) = uint8(full(X*(programBank(programNumber,:,programPointer))'));
            end

         

        end


        registerPointer = registerPointer + 1;





        if programSequenceLength ~= 1
            programPointer = programPointer + 1;
        end

    end
end