function [T] = thawTensor(fT,n)
    try
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % init return
        % there can be multiple columns of identically
        % frozen tensors - this will thaw all
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        T = {};
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get pointer list
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [ptrList,dm,type] = calculatePTRs(fT(:,1));
        % if nargin is one, then thaw all
        if nargin == 1
            n = 1:size(ptrList,2);
        end
        % for each grade
        for g = 1:numel(n)
            % get the first tensor for preallocation
            %V = fT(ptrList(1,n(g)):ptrList(2,n(g)),1);
            % for each tensor
            V = fT(ptrList(1,n(g)):ptrList(2,n(g)),:);
            % restore size
            V = reshape(V,[dm{n(g)}' size(V,2)]);

            % handle string case
            if (type{g} == 1)
                for e = 1:size(V,2)
                    tV{e} = char(V(:,e)');
                    tV{e} = strtrim(tV{e});
                end
                if numel(tV) == 1
                    V = tV{1};
                else
                    V = tV;
                end
            end



            % store output
            T{g} = V;
        end
        if numel(n) == 1
            T = T{1};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    catch ME
        ME
    end
    
end

function [ptrList,dm,type] = calculatePTRs(fT)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fT := folded tensor
    % ptrList := array point pointers [[str;stp]...[]]
    % dm := dimensions 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % measure grade
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ptrList = [];
    % grade is first in the order and the
    % the number of cells in the array
    grade = fT(1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for each rank
    ptr = 2;
    for g = 1:grade
        % get rank
        rnk = fT(ptr);
        % get dim of tensor
        dm{g} = fT(ptr+1:ptr+rnk);
        % get the type
        type{g} = fT(ptr+rnk+1);
        % set str -- added one more for type
        str = ptr+rnk+1+1;
        % set stp
        stp = ptr+rnk+1+prod(dm{g});
        % stack pointer list
        ptrList = [ptrList [str;stp]];
        % increment
        ptr = stp + 1;
    end
end
%{

% stack tensor to be frozen - cell format = direct sum
T = {rand(4,5),rand(3,4,5,6)};

fT = freezeTensor(T);
rT = thawTensor(fT);
rT{1} == T{1}
rT{2} == T{2}


T1 = {rand(4,5),rand(3,4,5,6),rand(1,2,3)};
T2 = {rand(4,5),rand(3,4,5,6),rand(1,2,3)};
fT1 = freezeTensor(T1);
fT2 = freezeTensor(T2);
singleOne = thawTensor([fT1]);
rTW = thawTensor([fT1;fT2]);
rTW{1} == T1{1}
rTW{2} == T1{2}



T1 = {{'a','b','c','d','e'},rand(4,5),rand(3,4,5),rand(1,2,5)};
fT1 = freezeTensor(T1,true);
singleOne = thawTensor([fT1]);



T1 = {rand(4,5),rand(3,4,5),rand(1,2,5)};
fT1 = freezeTensor(T1,true);
rTW = thawTensor(fT1);



%}













