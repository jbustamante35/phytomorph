function [I,GR] = groundRemove(I,GR)

    for j = 1:size(I,2)
        if nargin == 1
            tmp = I(:,j) > 0;
            tmp = [ones(size(tmp)) tmp ones(size(tmp))];
            tmp = imfill(logical(tmp),[1,2],8) - tmp;
            tmp = tmp(:,2);
            sidx = find(tmp);
        else
            sidx = -GR(j);
        end
        
        if ~isempty(sidx)
            I(:,j) = circshift(I(:,j),-sidx(end));
            GR(j) = -sidx(end);
        end
       
    end

end