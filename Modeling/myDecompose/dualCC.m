function [d] = dualCC(d,z)

        % get creator and carver
        w = z(:,1:size(d,2));
        v = z(:,(size(d,2)+1):end);

        
        d = createOP(d,w);
        d = carverOP(d,v);

end