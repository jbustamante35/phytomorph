function [G,img] = countOP2(data,vec,N)

    gmm = gmdistribution(vec(1:2)',reshape(vec(3:4)',[1 1 2]));
    vec = vec(5:end);
    
    
    for e = 1:size(data,3)
        source = mtimesx(vec,data(:,:,e));
        %source = vec'*data(:,:,e);
        %source = logsig(source);
        img = source;
        
        [P,nlogl] = posterior(gmm,source);
        
        
        
        
        %{
        if e == 1
            vS = im2col(source,[51 51],SZ);
        end
        %}
        
        %source = sort(source);
        %G(e) = norm(target - source');
        
        %CNT(e) = sum((source > .5).*source); 
        
    end
    
    G = norm(G);
    G2 = norm(CNT - N*ones(size(CNT)));
    %fprintf([num2str(G) '--' num2str(G2) '\n']);
    %mag
    G = G + G2;
end