function [STORE,DIS] = findProgram2(dispV,ML,LOOP,Y,type,X,MUTATE_NUMBER,MAX_WID,BYTE_SIZE,MAX_P_NUMBER,PROGRAM)
    %%%%%%%%%%%%%%%%%%%%%%%%
    % disp :- to display the results as program search proceeds
    % ML :- number of start points - loop over
    % LOOP :- iterations per start point
    % Y :- output variable looking to simuate
    % X :- inputs
    % MUTATE_NUMBER :- 

    evalProgram = @(x)evalSimMetrics(x,Y,type);
    
    % if disp
    if ~isempty(dispV)
        h1 = dispV(1);
        h2 = dispV(2);
        disp = true;
    else
        disp = false;
    end
    
    eP = [];
    delta = [];
    %%%%%%%%%%%%%%%%%%%%%%%%
    % loop over the number of start-points
    for ml = 1:ML
        
        %{
        vars = 4^2;
        degreeFreedom = nchoosek(16,0);
        degreeFreedom = nchoosek(8,1);
        degreeFreedom = nchoosek(4,2);
        degreeFreedom = nchoosek(2,3);
        degreeFreedom = nchoosek(1,4);
        degreeFreedom = nchoosek(0,5);
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%
        % randomly select program for init
        % need to change this to accommodate 2 bits
        % 255 needs to change
        %%%%%%%%%%%%%%%%%%%%%%%%
        PROGRAM_TEST = randi(MAX_P_NUMBER,1,MAX_WID,'uint8');
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % number of iterations per start program
        for loop = 1:LOOP
            %%%%%%%%%%%%%%%%%%%%%%%%
            % if the 'real' program is given and disp is on
            if nargin > 9
                delta(ml,loop) = bitwise_hamming(PROGRAM_TEST',PROGRAM');
                if disp
                    figure(h1);
                    plot(delta');
                end
            else
                delta(ml,loop) = 0;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%


            
            toMut = PROGRAM_TEST;
            cnt = 1;
           
            %%%%%%%%%%%%%%%%%%%%%%%%
            % init program mutation - get mutation numbers for random mutation(s)
            % number of bytes to mutate
            NUM_BYTE_TO_MUT = min(MUTATE_NUMBER,size(PROGRAM_TEST,2));
            % get the bytes to mutate at random
            rand_MUT_byte = randi(size(toMut,2),NUM_BYTE_TO_MUT,1);
           
            
            % for stochastic gradient decent - mut each bit in rand byte - type I
            bitMUTNumber = min((2^BYTE_SIZE),MAX_P_NUMBER);

            %{
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
            %%%%%%%%%%%%%%%%%%%%%%%%
            %}
            
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
            



            %{
            % for timing
            %fprintf(['Starting eval of programs @ ' num2str(size(progB,1)) '.\n']);
            tm = clock;
            preY = bo(X,progB(1,:));
            ftm = etime(clock,tm);
            %fprintf(['Estimated delta-T @ ' num2str(size(progB,1)*ftm/60) ' min.\n']);
            tm = clock;
            %}

            %%%%%%%%%%%%%%%%%%%%%%%%
            % start program evaluation
            % eval the mutated program(s) over the inputs
            preY = boo(X,progB);
            % evaluate the metric for the result of each program
            ME = zeros(size(progB,1),1);
            % loop over each program
            for iter = 1:size(progB,1)
                ME(iter) = evalProgram(preY(:,iter));
                %ME(iter) = -mutualinfo(preY(:,iter),Y);
                %ME(iter) = bitwise_hamming(preY(:,iter),Y);
            end
            % end program evaluation
            %%%%%%%%%%%%%%%%%%%%%%%%
            
            

            %{
            ME = zeros(size(progB,1),1);
            for iter = 1:size(progB,1)
                preY = bo(X,progB(iter,:));
                if ~all(preY == preY2(:,iter))
                    break
                end
                ME(iter) = mutualinfo(uint8(preY),Y);
            end
            %}

            ftm = etime(clock,tm);
            %fprintf(['Starting eval of programs @' num2str(size(progB,1)) '.\n']);
            %fprintf(['Real delta-T @ ' num2str(ftm/60) ' min.\n']);
            

            %initY = bo(X,toMut);
            %initE = mutualinfo(uint8(initY),Y);


            %[finalE,mi] = max(ME);
            % select the new program from the mutated
            [finalE,mi] = min(ME);
            
            
            % eval the metric between the new program and the output
            %finalE = bitwise_hamming(boo(X,PROGRAM_TEST),Y);
            %finalE = evalProgram(boo(X,progB(mi,:)));
            
            
            if loop > 1
                if finalE < eP(ml,loop-1)
                    % set focus to the new mutated program
                    PROGRAM_TEST = progB(mi,:);
                end
            else
                PROGRAM_TEST = progB(mi,:);
            end
            
            finalE = evalProgram(boo(X,PROGRAM_TEST));
            
            eP(ml,loop) = finalE;

            % display if needed
            if disp
                figure(h2);
                plot(eP');
                drawnow
            end
        end
        STORE(ml,:) = PROGRAM_TEST;
        DIS(ml,1) = eP(ml,end);
        DIS(ml,2) = delta(ml,end);
    end
    if disp
        %close(h1)
        %close(h2)
    end
end