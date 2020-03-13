function [G,bd,CNT,img] = countOP3(data,vec,N,T,P,sz)

    try
        if nargout > 1
            bd = [];
        end

        W = reshape(vec(1:prod(sz)),sz);
        vec(1:prod(sz)) = [];
        b = vec(1:sz(1));
        vec(1:sz(1)) = [];
        mag = vec(1:sz(1));
        vec = [];


        % normalize the projections and mag them
        for e = 1:size(W,1)
            W(e,:) = W(e,:) / norm(W(e,:));
            WW(e,:) = W(e,:);
            %W(e,:) = W(e,:) * mag(e);
        end
        % obtain angles
        L = tril(WW*WW');
        L = (sum(abs(L(:))) - size(W,1));
        
        
        % return output
        if nargout == 2
            source = mtimesx(W,data(:,:,1));
            %source = bsxfun(@plus,b',source);
            source = exp(source);
            source = bsxfun(@times,sum(source,1).^-1,source);
            img = source;
        end


        for e = 1:size(data,3)
            
            
            % make the distributions
            source = mtimesx(W,data(:,:,e));
            
            %{
            %source = bsxfun(@plus,b',source);
            source = exp(source);
            source = bsxfun(@times,sum(source,1).^-1,source);
            %}
            
            for i = 1:size(source,1)
                source(i,:) = normpdf(source(i,:),b(i),mag(i));
            end
            source = bsxfun(@times,sum(source,1).^-1,source);
            
            
            
            
            cnt = ((1:size(source,1))-1)*source;


            if ~isempty(T)
                %tGOAL(e) = norm(cnt - T(e,:))/numel(cnt);
                g = abs(cnt - T(e,:));
                g = sort(g,'descend');
                tGOAL(e) = mean(g(1:100));
                %tGOAL(e) = mean(g);
                %tGOAL(e) = max(abs(cnt - T(e,:)));
            else
                tGOAL(e) = 0;
            end


            %CNT(e) = sum(((1:size(source,1))-1)*source);
            CNT(e) = sum(cnt);

            %CNT(e) = sum((source > .5).*source); 

        end



        CNT(isnan(CNT)) = Inf;

        WTS = [1 10 100];
         WTS = [0 0 100];
        %WTS = [1 0 100];
        G = mean(abs(CNT-N));
        bd = [G L mean(tGOAL)];
        

        G;
        G = WTS*bd';
        bd = bd.*WTS;
        if G == 25
            stop = 1;
        end
        G;
        %G = norm(G);
        %G2 = norm(CNT - N*ones(size(CNT)));
        %fprintf([num2str(G) '--' num2str(G2) '\n']);
        %mag
        %G = G + G2;
    catch ME
        ME = 0;
    end
end