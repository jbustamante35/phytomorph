classdef dataPointSelector < handle
    
    
    properties
        dataPoints = {};
        dataAxisNames = {};
        
        
        %% total draws
        totalN = 100;
    end
    
    methods
    
        function [obj] = dataPointSelector()


        end
        
        function [] = setN(obj,N)
            obj.totalN = N;
        end


        function [] = addAxis(obj,X,name)
            obj.dataPoints{end+1} = X;
            obj.dataAxisNames{end+1} = name;
        end


        function [idx] = rDrawN(obj,N)
            for e = 1:numel(obj.dataPoints)
                sz(e) = numel(obj.dataPoints{e});
            end
            idx = randi(prod(sz),N,1);
            obj.totalN = obj.totalN - N;
        end


        function [idx] = rDrawN_fix1(obj,N,M)
            idx = [];

            for e = 1:numel(obj.dataPoints)
                sz(e) = numel(obj.dataPoints{e});
            end


            fidx1  = randi(prod(sz(1)),N,1);
            
            for e = 1:N
                fidx2 = randi(prod(sz(2:end)),M,1);
                lidx = sub2ind([sz(1) sz(2:end)],repmat(fidx1(e),[M 1]),fidx2);
                idx = [idx;lidx];
            end
    

        end


        function [dr] = drawToValue(obj,idx)
            for e = 1:numel(obj.dataPoints)
                sz(e) = numel(obj.dataPoints{e});
            end
            for e = 1:numel(sz)
                X{e} = cell(1);
            end
            [X{:}] = ind2sub(sz,idx);
            dr = {};
            for r = 1:numel(idx)
                for e = 1:numel(sz)
                    dr{r}{e} = obj.dataPoints{e}{X{e}(r)};
                end
            end
        end
    end
    
    
    methods (Static)
        function [out] = mat2cell(in)
            for e = 1:size(in,1)
                out{e} = in(e,:);
            end
        end
    end
end