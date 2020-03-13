%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxB = 8;
programSize = randi(2^maxB-1,1);
PROGRAM_BITS = randi(2^maxB-1,1,programSize,'uint32');
PROGRAM_BITS = unique(PROGRAM_BITS');

numberPureStates = randi(100,1,'uint32');
numberPureStates = 100;
totalStates = [];
for state = 1:numberPureStates
    randState = randi(2^maxB-1,1);
    newState = [randState,state];
    totalStates = [totalStates;newState];
end
X = qBit(totalStates);
PROGRAM = qBit([ones(size(PROGRAM_BITS)) PROGRAM_BITS]);
Y = PROGRAM*X;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate program bank
bitBank = [];
numberPrograms = 20;
for b = 1:numberPrograms
    programSize(b) = randi(2^maxB-1,1);
    PROGRAM_BITS = randi(2^maxB-1,programSize(b),1,'uint32');
    %bitBank = [bitBank;[ones(size(PROGRAM_BITS)) b*ones(size(PROGRAM_BITS)) PROGRAM_BITS]];
    bitBank = [bitBank;[b*ones(size(PROGRAM_BITS)) PROGRAM_BITS]];
end
bitBank = unique(bitBank,'rows');
programBank = qBit(bitBank);
%%
his = [];
MM = [];
for loop = 1:500
   
    
    yPre = programBank*X;

    
    
    
    score = [];
    numberPrograms = programBank.data;
    numberPrograms = max(numberPrograms(:,1));
    checkTime = [];
    yData = Y.data;
    pData = yPre.data;
    for program = 1:numberPrograms
        
        tic;
        %{
        delta = [];
        for e = 1:numberPureStates
            delta(e) = Y(1,e) == yPre(program,e);
        end
        %}
        tmpY = pData(pData(:,1) == program,2);
        delta2 = numberPureStates - numel(setdiff(yData(:,2),tmpY)) - numel(setdiff(tmpY,yData(:,2)));
        
        
        %sum(delta) == delta2;
        
        %score(program) = sum(delta);
        score(program) = delta2;
        checkTime(program) = toc;
    end
    fprintf(['Checked ' num2str(numberPrograms) ' programs @' num2str(mean(checkTime)) ':' ...
        num2str(mean(checkTime)*double(numberPrograms)) '\n']);
    
    
    [his(loop),sidx] = max(score);
    
    if max(his) <= his(loop)
        bestProgram = programBank(sidx,:);
        size(bestProgram.data,1);
    end
    programBank = programBank(sidx,:);

    %{
    q1 = bestProgram.data;
    q2 = PROGRAM.data;
    MM(loop) = numel(intersect(q1(:,2),q2(:,2)));
    %}
    
    %randDN = randi(size(programBank.data,1),1,'uint32');
    tic;
    randDN = 40;
    randDN = randi(size(programBank.data,1),randDN,1,'uint32');
    fprintf(['Generated ' num2str(size(randDN,1)) ' deletions in :' num2str(toc) '\n'])

    toAdd = 40;
    randAN = randi(2^maxB-1,toAdd,1);
    %randAN = FIX(1,2);
    
    
    % get the prorgram
    guts = programBank.data;
    
    guts(:,1) = 1;
    Orginal_guts = guts;
    
    
    Orginal_guts;
    
    for g = 1:numel(randDN)
        newProgram = Orginal_guts;
        newProgram(randDN(g),:) = [];
        newProgram(:,1) = 1 + g;
        guts = [guts;newProgram];
    end
    
    for g = 1:numel(randAN)
        newProgram = Orginal_guts;
        newProgram(:,1) = numel(randDN) + 1 + g;
        newProgram = [newProgram;[numel(randDN) + 1 + g randAN(g)]];
        
        guts = [guts;newProgram];
    end
    
    loop;
    %{
    guts(randDN,:) = [];
    guts = [guts;[repmat(guts(1,1),[numel(randAN) 1]) randAN]];
    %}
    
    guts = unique(guts,'rows');
    numberPrograms = 111;
    
    
    
    programBank = qBit(guts);
    plot(his)
    drawnow
end


%%




objA = qBit([A,ones(size(A))]);
B = randi(2^32-1,400,1,'uint32');
B(234) = A(1);
B(235) = A(1);
objB = qBit([ones(size(B)),B]);
objC = objB*objA;