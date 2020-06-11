%classdef basisT <  handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
classdef basisT <  doid & matlab.mixin.Copyable

    properties
        % "point"
        P;
        % basis
        E;
        % size
        Z;
        noCache = true;
    end
    
    methods
        
        % constructor for basis tensors
        function [obj] = basisT(E,P,Z)
            if nargin == 0;P = [];Z = [];E = [];end
            if nargin == 1;P = 0;Z = [];end
            obj.P = P;
            obj.E = E;
            obj.Z = Z;
        end
        
        function [] = apply(this,op,level)
            
        end
        
        function [] = applyDimOperation(this,op)
            for e = 1:size(this,1)
                a = a + this;
            end
        end
        
        function [C] = applyOperation(this,op,operand)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % unary operator
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nargin == 2
                % apply op over each of this
                for e = 1:numel(this)
                    this(e).E = op(this(e).E);
                end
                C = this;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % binary operator
            % three types:
            % 1) pairwise
            % 2) combination
            % 3) bsxfun
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif nargin == 3
                
                % init pairwise to false
                pairwise2 = false;
                % two types of subtract
                sz{1} = size(this);
                sz{2} = size(operand);
                % if number of elements in are equal
                pairwise1 = numel(sz{1}) == numel(sz{2});
                % if equal then 
                if pairwise1;pairwise2 = all(sz{1} == sz{2});end
                
                pairwisezip = pairwise1 && pairwise2;
                
                if pairwisezip
                    for a = 1:size(numel(this))
                        data = op(this(a).E,operand(a).E);
                        C(a) = basisT(data,this(a).P,this(a).Z);
                    end
                    C = reshape(C,size(this));
                else
                    % compute the size of the object
                    [m,largei] = max([numel(sz{1}) numel(sz{2})]);
                    smalli = setdiff(1:2,largei);
                    pd = ones(1,numel(sz{largei}) - numel(sz{smalli}));
                    sz{smalli} = [sz{smalli} pd];
                    newSZ = max([sz{smalli};sz{largei}],[],1);
                    
                    
                    for a = 1:numel(this)
                        for b = 1:numel(operand)
                            % compute data
                            data = op(this(a).E,operand(b).E);
                            % use the (P,Z) from A
                            C(a,b) = basisT(data,this(a).P,this(a).Z);
                        end
                    end
                    
                    
                    C = reshape(C,newSZ);
                end
            end
        end
        
        
        %{
        function [] = permute(this,order)
            pfunc = @(x)permute(x,order);
            this.applyOperation(pfunc);
        end
        %}

        %{
        function [] = reshape(this,order)
            pfunc = @(x)permute(x,order);
            this.applyOperation(pfunc);
        end
        %}


        function [C]  = transpose(this)
            C = copy(this);
            C.E = C.E.';
        end
        
        function [C]  = ctranspose(this)
            C = copy(this);
            C.E = C.E';
        end
        
        function [C] = mtimes(A,B)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if A is numeric - then using tensor notation to multiply the
            % data. this is done in the test script.
            % A(i,j,k)*B(l,k)=C(i,j,l) (sum over k);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isnumeric(A)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % get the X and Y data
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                X = B(1).E;
                Y = B(2).E;
                szX = size(X);
                szY = size(Y);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % bring the index to the front
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                nX = [A(1,:) setdiff(1:ndims(X),A(1,:))];
                nY = [A(2,:) setdiff(1:ndims(Y),A(2,:))];
                X = permute(X,nX);
                Y = permute(Y,nY);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % get the size pre-reshape
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                szX2 = size(X);
                szY2 = size(Y);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % reshape
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                tmpZ = size(A,2);
                X = reshape(X,[prod(szX2(1:tmpZ)) prod(szX2((tmpZ+1):end))]);
                Y = reshape(Y,[prod(szY2(1:tmpZ)) prod(szY2((tmpZ+1):end))]);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % multiply
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Z = mtimesx(X,'T',Y);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % restore shape
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Z = reshape(Z,[szX2((tmpZ+1):end) szY2((tmpZ+1):end)]);
                % construct return
                C = basisT(Z);
            else
                C = A.applyOperation(@(X,Y)mtimesx(X,Y),B);
                C.Z = B.Z;
            end
           
        end
        
        function [C] = plus(A,B)
            C = A.applyOperation(@(X,Y)plus(X,Y),B);
        end
        
        function [C] = minus(A,B)
            C = A.applyOperation(@(X,Y)minus(X,Y),B);
        end

        function [res] = norm(A)
            nfunc = @(x)norm(x(:));
            C = A.applyOperation(nfunc);
            res = C.E;
        end

        function [ret] = eq(A,B)
            C = A.applyOperation(@(X,Y)eq(X,Y),B);
            szC = size(C);
            for e = 1:numel(C)
                ret(e) = all(C(e).E(:));
            end
        end






        
        function [C] = or(this,A)
            C = this*A;
        end
        
        function [h] = gt(this,A)
            here = 1;
        end
        
        function [C] = sum(this,dim)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % currently this makes a copy of the data and operates on the
            % copy
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % bring the index to the front
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            nX = [dim setdiff(1:ndims(this),dim)];
            X = permute(this,nX);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get the size pre-reshape
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            szX = size(X);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % reshape
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tmpZ = 1;
            X = reshape(X,[prod(szX(1:tmpZ)) prod(szX((tmpZ+1):end))]);
            completeSZ1 = [1 szX((tmpZ+1):end)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % operation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for c = 1:size(X,2)
                id = 0; % the id of +
                for e = 1:size(X,1)
                    id = id + X(e,c).E;
                end
                C(c) = basisT(id);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            C = squeeze(reshape(C,[completeSZ1]));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % function level operations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        %{
        % MOVED TO BASISF
        % interpolate 2d-function
        % multi-purpose as cropping
        % x isa domain/pointset xor a box
        function [res] = f(this,x,layer)
            % copy to not over write the input
            res = copy(this);

            % switch on size for crop box
            if all(size(x(1).E) == [1 4])
                for e = 1:numel(res)
                    res(e).E = imcrop(res(e).E,x.E);
                end

            else
            

                if nargin == 2;layer = 1:size(res(1).E,3);end
                % init the size of the 
                sz = x.Z;
                % for each object in this
                for e = 1:numel(res)
                    % if this is empty or no-cache is true
                    if isempty(res(e).E) || (res.noCache)
                        % preallocate
                        if ~isempty(sz);patchSZ = [sz(1:2) numel(layer)];end
                        % interpolate
                        res(e).E = squeeze(ba_interp2(res.E,x.E(1,:),x.E(2,:)));
                        % reshape
                        if ~isempty(sz);res(e).E = reshape(res(e).E,patchSZ);end
                    end
                end
            end
        end
        %}


        % show as image-function
        function [] = imshow(this,idx)
            if nargin == 1;idx = 1:size(this,3);end
            imshow(this.E(:,:,idx),[]);
            drawnow
        end

        function [ret] = point2box(this,boxDim)
            ret = copy(this);
            for e = 1:numel(this)
                ret(e).E = point2Box(ret(e).E(1:2)',boxDim);
            end
        end


        %{
       
        
        function [E] = distance(A,B)
            for a = 1:numel(A)
                for b = 1:numel(B)
                    E(a,b,1) = A(a).displacement(B(b));
                    E(a,b,2) = A(a).stretch(B(b));
                    E(a,b,3) = A(a).rotation(B(b));
                end
            end
            E = squeeze(E);
        end
        
        function [totalD,displaceM] = displacement(A,B)
            totalD = norm(B.T(:,end) - A.T(:,end));
            if nargout == 2
                displaceM = eye(size(A.T,1));
                for e = 1:size(A.T,1)
                    displaceM(e,end) = B.T(e,end) -  A.T(e,end);
                end
                displaceM(end,end) = 1;
            end
        end
        
        function [totalS,stretchM] = stretch(A,B)
            for e = 1:(size(A.T,2)-1)
                %S(e) = norm(B.T(:,e)) - norm(A.T(:,e));
                S(e) = norm(B.T(:,e))/norm(A.T(:,e)) - 1;
            end
            % if asked return the stretch
            if nargout == 2
                stretchM = eye(size(A.T,1));
                for e = 1:numel(S)
                    stretchM(e,e) = S(e)+1;
                end
            end
            % return thte total stretch
            totalS = norm(S);
        end
        
        function [R] = rotation(A,B)
            for e = 1:(size(A.T,2)-1)
                vecA = A.T(:,e);
                vecB = B.T(:,e);
                vecA = vecA / norm(vecA);
                vecB = vecB / norm(vecB);
                R(e)= acos(vecA'*vecB);
            end
            R = norm(R);
        end
        
        %}
        
        
    
    end

    methods (Sealed)
    
        function [res] = prod(this,dim)
            res = copy(this(end));
            for e = 1:(numel(this)-1)
                res = this(end-e)*res;
            end
        end

    end
    
    methods (Static)
    
        function [A] = affine2D(P,T)
            toTranspose = false;
            if nargin == 1;T = eye(numel(P)+1);end
            %if nargin == 1;T = eye(numel(P)+1);end
            %if nargin == 1;toTranspose = false;end
            P = [P;1];
            %T = cat(2,T,P);
            T(:,end) = P;
            if toTranspose;T = T';end
            A = basisT(T,P,[]);
        end
        
        function [domain] = affineSpace2D(X,Y)
            [d1,d2] = ndgrid(X,Y);
            x = [d2(:),d1(:),ones(size(d1(:)))]';
            p = [0;0;1];
            domain = basisT(x,p,size(d1));
        end
        
        function [domain] = affineSpace2D_disk(T,R)
            [t,rad] = ndgrid(T,R);
            d1 = rad.*cos(t);
            d2 = rad.*sin(t);
            x = [d2(:),d1(:),ones(size(d1(:)))]';
            p = [0;0;1];
            domain = basisT(x,p,size(d1));
        end
        
        function [T] = transform2d2(varargin)
            if nargin == 3
                V(1) = varargin{1}(1);
                V(2) = varargin{1}(2);
                V(3) = varargin{2}(1);
                V(4) = varargin{3}(1);
                V(5) = varargin{3}(2);
                V(6) = varargin{3}(2);
            else
                V = cell2mat(varargin);
            end
            % displacement
            T1 = eye(3);
            T1(1,end) = V(1);
            T1(2,end) = V(2);
            % rotate
            T2 = [[cos(V(3)) sin(V(3)) 0];[-sin(V(3)) cos(V(3)) 0];[0 0 1]];
            T4 = [[cos(V(6)+V(3)) -sin(V(6)+V(3)) 0];[sin(V(6)+V(3)) cos(V(6)+V(3)) 0];[0 0 1]];
            % stretch
            T3 = zeros(3);
            T3(end,end) = 1;
            T3(1,1) = V(4);
            T3(2,2) = V(5);
            %T = T1*T2*T3;
            T = T1*T4*T3*T2;
            T = basisT(T,T(:,end),[]);
        end
        
        function [T] = transform2d(varargin)
            if nargin == 3
                V(1) = varargin{1}(1);
                V(2) = varargin{1}(2);
                V(3) = varargin{2}(1);
                V(4) = varargin{3}(1);
                V(5) = varargin{3}(2);
            else
                V = cell2mat(varargin);
            end
            % displacement
            T1 = eye(3);
            T1(1,end) = V(1);
            T1(2,end) = V(2);
            % rotate
            T2 = [[cos(V(3)) -sin(V(3)) 0];[sin(V(3)) cos(V(3)) 0];[0 0 1]];
            % stretch
            T3 = zeros(3);
            T3(end,end) = 1;
            T3(1,1) = V(4);
            T3(2,2) = V(5);
            %T = T1*T2*T3;
            T = T1*T2*T3*inv(T2);
            T = basisT(T,T(:,end),[]);
        end

        function [func,initGuess,mf] = makeTF(dX_range,dR_range,dS_range,initGuess)
            if nargin == 3;initGuess = [0 0 0 1 1];end
            initGuess(1) = (initGuess(1) - dX_range(2,1))*(dX_range(1,1) - dX_range(2,1))^-1;
            initGuess(2) = (initGuess(2) - dX_range(2,2))*(dX_range(1,2) - dX_range(2,2))^-1;
            initGuess(3) = initGuess(3);
            initGuess(4) = (initGuess(4) - dS_range(2,1))*(dS_range(1,1) - dS_range(2,1))^-1;
            initGuess(5) = (initGuess(5) - dS_range(2,2))*(dS_range(1,2) - dS_range(2,2))^-1;
            
            glf = @(X,M,m)((abs(X)>=1)*(sign(X)+(X==0))*M) + ...
                          (abs(X)<1)*(((sign(X)+(X==0))*(M-m)*abs(X))+(sign(X)+(X==0))*m);
                      
            glf = @(X,M,m)((abs(X)>=1)*(sign(X)+(X==0))*M) + ...
                          (abs(X)<1)*((M-m)*abs(X)+m);
           
            dX1 = @(X)glf(X,dX_range(1,1),dX_range(2,1));
            dX2 = @(X)glf(X,dX_range(1,2),dX_range(2,2));
            
            dR1 = @(X)glf(X,dR_range(1,2),dR_range(2,1));
            dR2 = @(X)glf(X,dR_range(1,2),dR_range(2,2));
            %dR = @(X)glf(X,dR_range(2),dR_range(1));
            %dR = @(X)pi*X;
            dS1 = @(X)glf(X,dS_range(1,1),dS_range(2,1));
            dS2 = @(X)glf(X,dS_range(1,2),dS_range(2,2));
            % stack functions
            F{1} = dX1;
            F{2} = dX2;
            F{3} = dR1;
            F{4} = dS1;
            F{5} = dS2;
            F{6} = dR2;
            % parameter mapping
            mf = @(X)arrayfun(@(f,X)F{f}(X),(1:numel(X)),X);
            % return generating function
            func = @(X)basisT.transform2d(mf(X));
        end
        


    end
    
    
    
end

%{
    clear e
    for i = 1:4
        for j = 1:9
            P = rand(1,3);
            E = rand(3,3);
            e(i,j) = basisT(P,[],E);
        end
    end

    % note: I is already loaded
    G = rgb2gray(I);
    [g1,g2] = gradient(G);
    globG = cat(4,G,g1,g2);
    fT = @(P,gD)transformFromGradient(gD,P);
    T = fwdT([30 30],fT);
    T.getT(globG);
   
    A = fwdT([40 50],fT);
    A.getT(globG);
    A.resetT();

    B = fwdT([40 51],fT);
    B.getT(globG);
    B.resetT();
    E = A.distance(B)


    B.T(1,1) = 2;
    E = A.distance(B)
    

    A.morph(B)
    clear pi
    alpha = pi/8;
    
    ROT = [[cos(alpha) sin(alpha) 0];[-sin(alpha) cos(alpha) 0];[0 0 1]];
    B.T = (A.T*ROT);
    morph(A,B)

%}