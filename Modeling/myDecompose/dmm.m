function [d,score] = dmm(z,d,u)
    try
        % turn z into logical
        z = round(z);
        z = bitget(z,1);
        
        % get creator and carver
        w = z(:,1:size(d,2));
        v = z(:,(size(d,2)+1):end);

        
        %u = createOP(u,w);
        %u = carverOP(u,v);
        
        
        %%%%%%%%%%%%
        % create only
        wu = createOP(u,w);
        d01 = pdist2(wu,d,'hamming');
        % carve only
        vu = carverOP(u,v);
        d10 = pdist2(vu,d,'hamming');
        % create then carve
        wuvu = carverOP(wu,v);
        d11 = pdist2(wuvu,d,'hamming');
        % no-thing
        d00 = pdist2(u,d,'hamming');
        d00 = repmat(d00,[size(z,1) 1]);
       
       
       
        [d,score] = min(cat(3,d00,d01,d10,d11),[],3);
        d = mean(d,2);
        
        
        s
    catch ME
        ME;
    end
end