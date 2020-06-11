function [d,score] = dm(z,d,u)
    try
        % turn z into logical
        z = round(z);
        z = bitget(z,1);
        % get creator and carver
        w = z(:,1:size(d,2));
        v = z(:,(size(d,2)+1):end);

        
        d = createOP(d,w);
        d = carverOP(d,v);
        
        
        %{
        % create only
        wu = bsxfun(@bitor,w,u);
        % carve only
        vu = bsxfun(@bitand,~v,u);
        % create then carve
        wuvu = bitand(wu,vu);
        % no-thing
        d00 = pdist2(u,d,'hamming');d00 = repmat(d00,[size(z,1) 1]);
        d01 = pdist2(wu,d,'hamming');
        d10 = pdist2(vu,d,'hamming');
        d11 = pdist2(wuvu,d,'hamming');
        [d,score] = min(cat(3,d00,d01,d10,d11),[],3);
        d = mean(d,2);
        %}
        
        s
    catch ME
        ME;
    end
end