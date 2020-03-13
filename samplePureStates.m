function [INT] = samplePureStates(MAX_B,TALL)
    % generate 'TALL' random numbers
    % randi generates from 1->2^MAXB
    % subtract 1 to be between 0->2^N-1
    % old
    %INT = double(randi((2^MAX_B-1),TALL,1)-1);
    % new
    INT = double(randi((2^MAX_B),TALL,1)-1); 
end