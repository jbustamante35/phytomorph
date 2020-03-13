function [X] = transformLength(X,,direction)
   
    X = X(:,1);

    if direction == -1
        l = bxsfun(@minus,X((:,end),l);
    else
        l = bxsfun(@minus,l(1,end),l);
    end
end