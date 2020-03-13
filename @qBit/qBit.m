classdef qBit < handle
    
    properties
        data;
    end
    
    methods
        function [obj] = qBit(data)
            obj.data = data;
        end
        
        function [] = permute(obj,p)
            if numel(p) == size(obj.data,2)
                obj.data = obj.data(:,p);
            else
                fprintf(['error in dims\n']);
            end
        end
        
        function [] = insertState(obj,newState)
            obj.data = [obj.data newState];
        end
        
        function [ret] = subsref(obj,subs)
            if strcmp(subs.type ,'()')
                %if numel(subs.subs) == size(obj.data,2)
                IDX = [];
                for i = 1:numel(subs.subs)
                    if strcmp(subs.subs{i},':')
                        IDX = [IDX ones(size(obj.data,1),1)];
                    else
                        IDX = [IDX any(bsxfun(@eq,obj.data(:,i),subs.subs{i}),2)];
                    end
                end
                %end
                if ~isempty(IDX)
                    IDX = all(IDX,2);
                end
                strip = obj.data(IDX,:);
                ret = qBit(strip);
            elseif strcmp(subs.type,'.')
                ret = builtin('subsref',obj,subs);
            end
        end
        
        function [ret] = eq(A,B)
            if isempty(A.data) && isempty(B.data)
                ret = 1;
            elseif numel(A.data) == numel(B.data)
                if  ~isempty(A.data) == ~isempty(B.data(:))
                     ret = 1;
                else
                     ret = 0;
                end
            else
                ret = 0;
            end
        end
        
        function [objC] = mtimes(objA,objB)
            toCompareA = size(objA.data,2);
            toCompareB = 1;
            dataC = [];
            for state = 1:size(objB.data,1)
                b = bsxfun(@eq,objA.data(:,toCompareA),objB.data(state,toCompareB));
                
                newCrow = objA.data(b,1:(end-1));
                %newCrow = unique(newCrow);
                newCcol = repmat(objB.data(state,2),[size(newCrow,1) 1]);
                dataC = [dataC;[newCrow newCcol]];
            end
            dataC = unique(dataC,'rows');
            
            objC = qBit(dataC);
        end
        
    end
end

%{

    A = randi(2^32-1,400,1,'uint32');
    objA = qBit([A,ones(size(A))]);
    B = randi(2^32-1,400,1,'uint32');
    B(234) = A(1);
    B(235) = A(1);
    objB = qBit([ones(size(B)),B]);
    objC = objB*objA;


%}