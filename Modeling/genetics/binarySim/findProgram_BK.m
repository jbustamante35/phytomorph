function [STORE,DIS] = findProgram(N,dispV,ML,LOOP,evalProgramOutputs,X,Y,MUTATE_NUMBER,MAX_WID,BYTE_SIZE,MAX_P_NUMBER,PROGRAM)
    %%%%%%%%%%%%%%%%%%%%%%%%
    % disp :- to display the results as program search proceeds
    % ML :- number of start points - loop over
    % LOOP :- iterations per start point
    % Y :- output variable looking to simuate
    % X :- inputs
    % MUTATE_NUMBER :- 

    
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
        
        [PROGRAM_TEST] = generateInitProgram(X,MAX_P_NUMBER,1,MAX_WID);
        
        %
        %{
        ridx = find(Y==1);
        PROGRAM_TEST = zeros(N,MAX_WID,'uint8');
        for iter = 1:numel(ridx)
            PROGRAM_TEST = bitor(PROGRAM_TEST,X(ridx(iter),:));
        end
        %}
        
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % number of iterations per start program
        for loop = 1:LOOP

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if the 'real' program is given and disp is on
            if nargin > 11

                delta(ml,loop) = measureProgramDistance(PROGRAM_TEST',PROGRAM');
                if disp
                    figure(h1);
                    plot(delta');
                end
            else
                delta(ml,loop) = 0;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            MUTATE_NUMBER = 10;
            [progB] = generateProgramBank(N,MUTATE_NUMBER,MAX_WID,BYTE_SIZE,MAX_P_NUMBER,PROGRAM_TEST,'full_gradient');            
            progB = progB{1};




            %%%%%%%%%%%%%%%%%%%%%%%%
            % start program evaluation
            % eval the mutated program(s) over the inputs
            [preY] = evalProgram(progB,X);
            % evaluate the metric for the result of each program
            ME = zeros(size(progB,1),1);
            % loop over each program
            for iter = 1:size(progB,1)
                if ~issparse(progB)
                    Ydata = preY(:,iter);
                else
                    Ydata = uint8(full(preY(:,iter)));
                end
                ME(iter) = evalProgramOutputs(Ydata);
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

            %ftm = etime(clock,tm);
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
            
            
            tmpY = evalProgram(PROGRAM_TEST,X);
            
            
            finalE = evalProgramOutputs(tmpY);
            
            eP(ml,loop) = finalE;

            % display if needed
            if disp
                figure(h2);
                plot((eP)');
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