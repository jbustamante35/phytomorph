function [G,img] = countOP(data,vec,N)
    mag = vec(end);
    vec(end) = [];
    b = vec(end);
    vec = vec(1:end-1);
    vec = vec / norm(vec);
    vec = vec * mag;
    target = zeros(size(data,2),1);
    target(1:N) = 1;
    target = sort(target);
    
    if nargout == 2
        source = mtimesx(vec,data(:,:,1)) + b;
        %source = vec'*data(:,:,e);
        source = logsig(source);
        img = source;
    end
    
    
    parfor e = 1:size(data,3)
        source = mtimesx(vec,data(:,:,e)) + b;
        %source = vec'*data(:,:,e);
        source = logsig(source);
        %img = source;
        
        %{
        if e == 1
            vS = im2col(source,[51 51],SZ);
        end
        %}
        
        source = sort(source);
        G(e) = norm(target - source');
        
        CNT(e) = sum((source > .5).*source); 
        
    end
    
    G = norm(G);
    G2 = norm(CNT - N*ones(size(CNT)));
    %fprintf([num2str(G) '--' num2str(G2) '\n']);
    %mag
    G = G + G2;
end