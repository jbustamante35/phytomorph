function [progBank] = generateProgramBank(N,MUTATE_NUMBER,MAX_WID,BYTE_SIZE,MAX_P_NUMBER,PROGRAM_TEST,generateType)
    for n = 1:size(PROGRAM_TEST,1)

            toMut = PROGRAM_TEST(N,:);
            cnt = 1;
           
            %%%%%%%%%%%%%%%%%%%%%%%%
            % init program mutation - get mutation numbers for random mutation(s)
            % number of bytes to mutate
            NUM_BYTE_TO_MUT = min(MUTATE_NUMBER,size(PROGRAM_TEST,2));
            % get the bytes to mutate at random
            rand_MUT_byte = randi(size(toMut,2),NUM_BYTE_TO_MUT,1);
            % for stochastic gradient decent - mut each bit in rand byte - type I
            bitMUTNumber = min((2^BYTE_SIZE),MAX_P_NUMBER);


            switch generateType
                case 'stochastic'
                    %%%%%%%%%%%%%%%%%%%%%%%%
                    % begin mutation
                    %fprintf(['Starting mutation of programs @ bytes @ ' num2str(NUM_BYTE_TO_MUT) '.\n']);
                    % create a new bank of mutated programs
                    progB = zeros(NUM_BYTE_TO_MUT*(2^BYTE_SIZE),size(PROGRAM_TEST,2),'uint8');
                    tm = clock;
                    for b = 1:NUM_BYTE_TO_MUT
                        byte = rand_MUT_byte(b);
                        for bit = 1:bitMUTNumber
                            newP = toMut;
                            newP(byte) = bitset(newP(byte),bit,not(bitget(newP(byte),bit)));
                            progB(cnt,:) = newP;
                            cnt = cnt + 1;
                        end
                    end
                    ftm = etime(clock,tm);
                    %fprintf(['Starting mutation of programs @ bytes @  ' num2str(ftm/60) ' min.\n']);
                     % end mutation
                     
                     
                case 'stochastic2'
                    %%%%%%%%%%%%%%%%%%%%%%%%
                    % begin mutation
                    %fprintf(['Starting mutation of programs @ bytes @ ' num2str(NUM_BYTE_TO_MUT) '.\n']);
                    % create a new bank of mutated programs
                    progB = zeros(NUM_BYTE_TO_MUT,size(PROGRAM_TEST,2),'uint8');
                    tm = clock;
                    for b = 1:NUM_BYTE_TO_MUT
                        byte = rand_MUT_byte(b);
                        newP = toMut;
                        
                        for bit = 1:bitMUTNumber
                            if randi(2) == 1
                                newP(byte) = bitset(newP(byte),bit,not(bitget(newP(byte),bit)));
                            end
                        end
                        
                        progB(cnt,:) = newP;
                        cnt = cnt + 1;
                    end
                    ftm = etime(clock,tm);
            
                case 'full_gradient'
                    
                    if ~issparse(PROGRAM_TEST)
                        % create a new bank of mutated programs
                        progB = zeros(NUM_BYTE_TO_MUT*bitMUTNumber,size(PROGRAM_TEST,2),'uint8');

                        % for gradient decent
                        tm = clock;
                        cnt = 1;
                        for byte = 1:MAX_WID
                            for bit = 1:bitMUTNumber
                                newP = toMut;
                                newP(byte) = bitset(newP(byte),bit,not(bitget(newP(byte),bit)));
                                progB(cnt,:) = newP;
                                cnt = cnt + 1;
                            end
                        end
                    else
                        % create a new bank of mutated programs
                        progB = sparse(size(PROGRAM_TEST,2),size(PROGRAM_TEST,2));
                        for e = 1:size(PROGRAM_TEST,2)
                            vi = ~logical(PROGRAM_TEST(1,e));
                            progB(e,:) = PROGRAM_TEST;
                            progB(e,e) = vi;
                        end
                    end
            end

        progBank{n} = progB;
    end
end