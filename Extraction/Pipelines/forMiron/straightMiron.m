function [funcY] = straightMiron(Sc1,dim,funcS,N)
    try
        disp = 0;
        emptyFlag = true;
        while emptyFlag
            %bins = linspace(min(Sc1(dim(1),:)),max(Sc1(dim(1),:)),N);
            [Y,bins] = discretize(Sc1(dim(1),:),N);
            emptyFlag = (numel(bins)-1) ~= numel(unique(Y));
            N = N - 1;
        end
        
        % check for empty bins
        

        
        for e = 1:(numel(bins)-1)
            fidx = find(Y == e);
            
            
            tmp = Sc1(:,fidx);
            
            
            for f = 1:numel(funcS)
                funcSY(f,e) = funcS{f}(tmp(dim(2),:));
            end
            
        end
        
        for f = 1:numel(funcS)
            yValues = funcSY(f,:);
            yValues = imfilter(yValues,fspecial('average',[1 11]),'replicate');
            funcY(f).func = @(values)interp1(bins(1:end-1),yValues,values,'linear','extrap');
            funcY(f).dimX = dim(1);
            funcY(f).dimY = dim(2);
        end
      
        
        
    catch ME
        ME
    end

end
