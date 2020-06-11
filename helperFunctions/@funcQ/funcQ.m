classdef funcQ < fwdT
    
    properties
        % sampled function at P
        F;
        globData;
        noCache = true;
    end
    
    methods

        function [obj] = funcQ(P,fT,globData)
            obj = obj@fwdT(P,fT);
            obj.globData = globData;
        end


        % sample the function
        % x - domain
        % sz - size of domain
        function [F] = getF(obj,x,layer,id)
            
            if nargin == 2;layer = 1:size(obj.globData,3);end
            
            if nargin < 4;id = 1;end
            sz = x.sz;

            F = [];
            for e = 1:numel(obj)
                if isempty(obj(e).F) || (obj.noCache)
                    if ~isempty(sz);patchSZ = [sz(1:2) numel(layer) numel(id)];end
                    T = obj(e).getT(obj.globData);
                    x.domain = mtimesx(T,x.domain,'T')';
                    v = zeros([size(x.domain,1) numel(layer) numel(id)]);
                    for i = 1:numel(id)
                        
                        v(:,:,i) = squeeze(ba_interp2(obj.globData(:,:,:,id(i)),x.domain(:,1),x.domain(:,2)));
                        
                        %{
                        for l = 1:numel(layer)
                            v(:,l,i) = ba_interp2(globData(:,:,layer(l),id(i)),x(:,1),x(:,2));
                        end
                        %}
                    end
                    obj(e).F = v;
                    if ~isempty(sz);obj(e).F = reshape(obj(e).F,patchSZ);end
                end
                if ~isempty(sz);F = cat(ndims(obj(e).F)+1,F,obj(e).F);end
            end
        end
        
        
        
        function [F] = F(obj,x)
        
            
        end
        
        
        function [] = distance(A,B)
            E = distance@fwdT(A,B);
            for a = 1:numel(A)
                for b = 1:numel(B)
                    A(a).F - B(b).F;
                end
            end
        end
        
        function [ret] = push(obj,dT)
            if nargout == 1
                ret = push@fwdT(obj,dT);
                ret.F = [];
             
            else
                push@fwdT(obj,dT);
                obj.F = [];
            end
        end
       
        function [d] = dY(A,B,x,t)
                a = A.getF(x);
                if nargin == 4;B.setT(reshape(t,size(A.T)));end
                b = B.getF(x);
                d = norm(a(:)-b(:));
        end

        function [tP,T] = morph(obj,target,x)
          
            %initF1 = obj.getF(x);
            %initF2 = target.getF(x);

            %solutionScale = 1000;
            %solutionOffset = 1000;
            %mapBack = @(t)(solutionScale^-1)*(t-solutionOffset);
            %mapForward = @(t)solutionScale*t + solutionOffset;

            initT = target.getT([]);
            func = @(t)obj.dY(target,x,t);

            %{
            swarmSize = 200;
            ops = optimoptions('particleswarm','UseParallel',true,'SwarmSize',swarmSize);
            per = .4;
            lb = [zeros(6,1);per*initT(7:8)';1];
            ub = [2*ones(6,1);(1+per)*initT(7:8)';1];
            T = particleswarm(func,numel(initT),lb,ub,ops);
            initT = T(:);
            %}
            %{
            per = .4;
            lb = [zeros(6,1);per*initT(7:8)';1];
            ub = [2*ones(6,1);(1+per)*initT(7:8)';1];
            T = ga(func,numel(initT),[],[],[],[],lb,ub)
            %}
            %{
            tic
            initT = T(:);
            ops = optimoptions('fminunc','Display','none','FiniteDifferenceType','central');
            T = fminunc(func,initT(:)',ops);
            toc
            %}
            per = .4;
            lb = [zeros(6,1);per*initT(7:8)';1];
            ub = [2*ones(6,1);(1+per)*initT(7:8)';1];

            ops = optimoptions('patternsearch','Display','none','UseParallel',true);
            T = patternsearch(func,initT(:)',[],[],[],[],lb,ub,[],ops);
            %tic
            %initT = T(:);
            %ops = optimset('Display','none');
            %T = fminsearch(func,initT(:)',ops);
            %toc

            T = reshape(T,size(obj.T));
            newP = T(1:2,3)';
            tP = funcQ(newP,obj.fT,target.globData);
            %target.T = finalT;
            %target.setP();
            %finalF2 = tP.getF(x);

        end
    

        function [tP,T] = morph2(obj,target,x)
            x1 = mtimesx(obj.T,x.domain,'T')';
            threshold = .001;
            I(:,:,1) = obj.getF(x);
            I(:,:,2) = target.getF(x);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % take the gradient
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [D1 D2] = gradient(I(:,:,1));
            dX = [1 1];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % interpolate
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            externalInterp = true;
            
            if externalInterp
                Ii = ba_interp2(I(:,:,1),x1(:,1),x1(:,2));
                D1i = ba_interp2(D1,x1(:,1),x1(:,2));
                D2i = ba_interp2(D2,x1(:,1),x1(:,2));
            else
                Ii = interp2(I(:,:,1),x1(:,1),x1(:,2));
                D1i = interp2(D1,x1(:,1),x1(:,2));
                D2i = interp2(D2,x1(:,1),x1(:,2));
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % interpolate
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            D = [D1i D2i];

            icnt = 1;
            flag = 1;
            N = [];

            init_T = target.T(1:2,:);
            init_T = init_T(:)';
            while flag && (norm(dX) > threshold)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % init transformation
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                TR = reshape(init_T,[2 3]);
                Xt = mtimesx(TR,x.domain,'T')';
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % if internal
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if externalInterp
                    Gi = ba_interp2(I(:,:,2),Xt(:,1),Xt(:,2));
                else
                    Gi = interp2(I(:,:,2),Xt(:,1),Xt(:,2));    
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % solution vector
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Mi = [D.*repmat(x.domain(:,1),[1 2]) D(:,1) D.*repmat(x.domain(:,2),[1 2]) D(:,2)];
                dY = Mi\(Ii-Gi);
                dX = [dY(1) dY(2) dY(4) dY(5) dY(3) dY(6)]';
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % displace vector
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                init_T = init_T + dX';
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % measure image distance
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                N(icnt) = norm(Ii(:)-Gi(:));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % ensure that the norm is minimizing.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if icnt >= 2
                    flag = N(icnt) <  N(icnt-1);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                icnt = icnt + 1;
            end
            T = [initT 0 0 1];
            T = reshape(T,[3 3]);
        end
    end
    
end

%{
    % note: I is already loaded
    G = rgb2gray(I);
    [g1,g2] = gradient(G);
    globG = cat(4,G,g1,g2);
    fT = @(P,gD)transformFromGradient(gD,P);
    T = fwdT([30 30],fT);
    T.getT(globG);



    [x1,x2] = ndgrid(linspace(-10,10,20),linspace(-10,10,20));
    szX = size(x1);
    x = [x1(:) x2(:) ones(size(x1(:)))];
    P = funcP([30 30],fT);
    F = P.getF(globG,x,szX);
    arrayP = [P;P];
    F = arrayP.getF(globG,x,szX);



    P.distance(copy(P));


%}