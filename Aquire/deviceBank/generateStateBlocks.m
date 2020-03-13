function [stateBlocks] = generateStateBlocks(sImage,portIdx,TotalNumber,SampleNumber,transitionFs,successFs,failFs)
    T(:,:,1) = eye(3);
    T(:,:,2) = [[0 1 0];[0 0 1];[1 0 0]];
    I = eye(2);
    TT = zeros(6);
    for e = 1:size(I,1)
         II = diag(I(e,:));
         TT = TT + kron(II,T(:,:,e));
    end
   
    
    % make matrix of event tiles
    stateBlocks = [];
    
    % make first tile in sequence 
    tF = transitionFs{1};
    sF = successFs{1};
    fF = failFs{1};
    curBlock = stateBlock(TT,tF,sF,fF,stateDirac(2,3));
    
    
    % make a state transistion filter
    lFilter = myStateTransitionListener([0;1;0],[1;0;0],@(X,Y)testListen);
    lFilter.addFilterListener(curBlock,'currentState');
 
    
    % make a state transistion filter
    lFilter = myStateTransitionListener([1;0;0],[0;0;1],@(X,Y)curBlock.setState([0;1;0]));
    lFilter.addFilterListener(curBlock,'currentState');
   
    
    
    dispFunc=@(msg)sImage.updateTile(portIdx,1,rgb,sat,subImage)
    
    
    stateBlocks = [stateBlocks,curBlock];
    for e = 2:(TotalNumber-1)
        
        if e <= SampleNumber + 1
            tF = transitionFs{2};
            sF = successFs{2};
            fF = failFs{2};
        end
        
        if e > SampleNumber
            tF = transitionFs{3};
            sF = successFs{3};
            fF = failFs{3};
        end
        
        stateBlocks = [stateBlocks,stateBlock(TT,tF,sF,fF,stateDirac(3,3))];
    end
    
    tF = transitionFs{4};
    sF = successFs{4};
    fF = failFs{4};
    stateBlocks = [stateBlocks,stateBlock(TT,tF,sF,fF,stateDirac(2,3))];
    
    
    
end