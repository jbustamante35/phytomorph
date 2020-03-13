function [ret,BEST] = sparseProgramFinder(X,Y,B,onesN,W,UC,disp)

    TOTV = std(Y).^2;
    
    
    PER = .5;
    testIDX = [];
    for e = 1:size(X,2)
        encoding(e) = find(X(:,e)==1);
    end
    UQ = unique(encoding);
    for u = 1:numel(UQ)
        sE(u) = sum(encoding==UQ(u));
        takeN = floor(sE(u)*PER);
        fidx = find(encoding == UQ(u));
        fidx = fidx(randperm(numel(fidx)));
        testIDX = [testIDX;fidx(1:takeN)'];
    end
    
    XBK = X;
    YBK = Y;
    
   %{
    per = round(PER*size(X,2));
  
    testIDX = randi(size(X,2),1,per);
    %testIDX = 1;
    %}
    
    
    trainIDX = setdiff(1:size(X,2),testIDX);
    
    testX = X(:,testIDX);
    X = X(:,trainIDX);
    
    
    testY = Y(testIDX);
    Y = Y(trainIDX);

    if disp
        h1 = figure;
        %h2 = figure;
    end


    % number of trials
    T = 1;
    % number of cycles
    CY = 1;
    % std of entropy
    STD = 1;
    % tot plot
    tot = [];
    % max
    mx = 0;
    for t = 1:T

        
        X = XBK;
        Y = YBK;
        
        %{
        per = round(PER*size(X,2));
        XBK = X;
        YBK = Y;
        testIDX = randi(size(X,2),1,per);
        %testIDX = 1;
        trainIDX = setdiff(1:size(X,2),testIDX);
        %}
        
        statesPerBin = size(X,1)/B;
        
        PER = .5;
        testIDX = [];
        for e = 1:size(X,2)
            encoding(e) = find(X(:,e)==1);
        end
        %UQ = unique(encoding);
        UQ = 1:size(X,1);
        for u = 1:numel(UQ)
            sE(u) = sum(encoding==UQ(u));
            takeN = floor(sE(u)*PER);
            fidx = find(encoding == UQ(u));
            fidx = fidx(randperm(numel(fidx)));
            testIDX = [testIDX;fidx(1:takeN)'];
            
            jump = setdiff(fidx,fidx(1:takeN));
            
            YVV(u) = mean(Y(jump));
            
        end
        %[~,yidx] = sort(YVV,'descend');
        [~,yidx] = sort(YVV,'ascend');
        
        
        trainIDX = setdiff(1:size(X,2),testIDX);
        
        
        testX = X(:,testIDX);
        X = X(:,trainIDX);


        testY = Y(testIDX);
        Y = Y(trainIDX);


        %%%%%%%%%%%%%%%%%%%%
        % init guess program
        %bs = randperm(B);
        bs = yidx;
        N = sparse(W,B);
        for prog = 1:(W+1)
            %N(prog,bs(1:onesN(prog))) = 1;
            N(prog,bs(1:statesPerBin)) = 1;
            bs(1:statesPerBin) = [];
        end


        % loop for gradient descent
        for cycle = 1:CY

            baseline = UC*N*X;
            baseline = norm(baseline - Y);
            rb = randi(size(N,2),100,1);
            %rb = 1:size(N,2);
            %%%%%%%%%%%%%%%%%%%%
            % for each prog
            for prog = 1:size(N,1)
               
                %
                %%%%%%%%%%%%%%%%%%%%
                % for each bit
                for b = 1:numel(rb)
                    bit = rb(b);
                    % init grad program
                    gN = N;
                    
                    bitValue = gN(prog,bit);
                    gN(:,bit) = 0;
                    gN(prog,bit) = ~bitValue;
                    
                    %gN(prog,bit) = ~gN(prog,bit);
                    %gY = (1:W)*gN*X;
                    gY = UC*gN*X;
                    %tmp = norm(gY - Y);
                    %tmp = -tmp/TOTV;
                    tmp = -corr(gY',Y');
                    if bitValue == 1
                        tmp = Inf;
                    end
                    %tmp = -std(gY - Y)^2 / TOTV;
                    %delta = normpdf(gY - Y,0,STD);
                    %grad(prog,bit) = -sum(log(delta));
                    grad(prog,b) = tmp;

                end
            end

            % change
            %{
            [V,sidx] = sort(grad(:),'ascend');
            
            NC = 1;
            for n = 1:NC
                [prog,b] = ind2sub(size(grad),sidx(n));
                bit = rb(b);
                if V(n) < baseline
                    bitValue = N(prog,bit);
                    N(:,bit) = 0;
                    N(prog,bit) = ~bitValue;
                end
            end
            %}
            %N(prog,bit) = ~N(prog,bit);



            gY = UC*N*X;
            if any(gY==0)
                here = 1;
            end
            %delta = norm(gY - Y);
            delta = -corr(gY',Y');
            %delta = std(gY - Y)^2 / TOTV;
            test_gY = UC*N*testX;
            %test_delta = norm(test_gY - testY);
            test_delta = -corr(test_gY',testY');
            %test_delta = std(test_gY - testY)^2 / TOTV;
            
            
            %delta = normpdf(gY - Y,0,STD);
            %tot(t,cycle) = -mean(log(delta));
            tot(t,cycle) = delta^2;
            tot(t,cycle+1:end) = max(tot(:));
            
            test_tot(t,cycle) = test_delta^2;
            test_tot(t,cycle+1:end) = max(test_tot(:));



            if tot(t,cycle) > mx
                mx = tot(t,cycle);
                gY = UC*N*X;
                
                
                if disp
                    figure(h1)
                    subplot(1,2,2)
                    plot(Y+.1*rand(size(Y)),gY+.1*rand(size(Y)),'b.')
                    hold on
                    plot(testY+.1*rand(size(testY)),test_gY+.1*rand(size(test_gY)),'ro')
                    hold off
                    %axis([0 8 0 8]);
                    drawnow


                    hold on
                  %  YPT = (1:W)*N*XT;
                  %  plot(YT,YPT,'ro')
                    hold off
                    drawnow
                end
                BEST = N;
            end

            
            
            if disp
                figure(h1)
                subplot(1,2,1);
                plot(tot')
                hold on
                plot(test_tot','r')
                hold off
                drawnow
            end


        end
    end
    %ret = mx;
    ret = [];
    ret(2) = mean(test_tot(:,end));
    ret(1) = mean(tot(:,end));
    if disp
        close(h1);
    end
    %close(h2);
end