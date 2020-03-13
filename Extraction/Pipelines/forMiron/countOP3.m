function [G,img] = countOP3(data,vec,N,sz)


    W = reshape(vec(1:prod(sz)),sz);
    vec(1:prod(sz)) = [];
    b = vec(1:sz(1));
    vec(1:sz(1)) = [];
    mag = vec(1:sz(1));
    vec = [];

    
    
    for e = 1:size(W,1)
        W(e,:) = W(e,:) / norm(W(e,:));
        WW(e,:) = W(e,:);
        W(e,:) = W(e,:) * mag(e);
    end
    
    L = tril(WW*WW');
    L = 10*(sum(abs(L(:))) - size(W,1));
    
    
    if nargout == 2
        source = bsxfun(@plus,b',mtimesx(W,data(:,:,1)));
        source = exp(source);
        source = bsxfun(@times,sum(source,1).^-1,source);
        img = source;
    end
    
    
    
    parfor e = 1:size(data,3)
        source = bsxfun(@plus,b',mtimesx(W,data(:,:,e)));
        source = exp(source);
        source = bsxfun(@times,sum(source,1).^-1,source);
        
        %CNT(e) = sum(((1:size(source,1))-1)*source);
        CNT(e) = sum((((1:size(source,1))-1)*source));
        
        %CNT(e) = sum((source > .5).*source); 
        
    end
    CNT(isnan(CNT)) = Inf;
    
    
    
    G = mean(abs(CNT-N));
    G;
    G = G + L;
    if G == 25
        stop = 1;
    end
    G;
    %G = norm(G);
    %G2 = norm(CNT - N*ones(size(CNT)));
    %fprintf([num2str(G) '--' num2str(G2) '\n']);
    %mag
    %G = G + G2;
end