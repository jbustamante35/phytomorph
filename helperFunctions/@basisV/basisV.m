classdef basisV <  handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable

    properties
        % "point"
        P;
        % basis
        E;
        % size
        Z;
        % basisGenerator
        eGen;
    end
    methods
        
        % constructor for basis tensors
        function [obj] = basisV(P,Z,eGen)
            if nargin == 0;Z = [];eGen = [];end
            obj.P = P;
            obj.E = [];
            obj.Z = Z;
            obj.eGen = eGen;
        end
        
        % compute the basis for each object in the array
        function [T] = computeE(obj,globData,force)
            if nargin == 2;force = false;end
            T = zeros(prod(obj(1).Z),numel(obj));
            % for each object in the array
            for e = 1:numel(obj)
                % 
                if isempty(obj(e).T) || force
                    obj(e).E = obj(e).eGen(obj(e).P,globData);
                end
                T(:,e) = obj(e).E;
            end
        end
        
        % hard set/force the T
        function [] = setE(obj,E)
            for e = 1:numel(obj)
                obj(e).T = E(:,:,e);
            end
        end
        
        
        % 
        function [] = setPT(obj,newP,globData)
            obj.P = newP;
            obj.getT(globData);
        end
        
        
        
        
        function [] = resetT(obj)
            szT = size(obj.T)-1;
            obj.T(1:szT(1),1:szT(2)) = eye(szT(1));
        end
        
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
        
        
        function [] = morph(A,B)
            [totalS,stretchM] = A.stretch(B);
            [totalD,displaceM] = A.displacement(B);
            Ap = A.T*stretchM*displaceM;
            %B.T = Ap*U;
            rotationM = B.T / Ap;
            target = Ap*rotationM;
        end
        
        
        function [ret] = push(A,dT)
            if nargout == 0
                %A.T = A.T*dT;
                A.T = mtimesx(A.T,dT);
            else
                ret = copy(A);
                %ret.T = ret.T*dT;
                ret.T = mtimesx(ret.T,dT);
            end
        end
        
        %{
        function [varargout] = subsref(obj,S)
             %varargout = cell(1,nargout);
             %[varargout{:}] = builtin('subsref',obj,S);
            
            
            if strcmp(S(1).type,'{}')
                if strcmp(S(1).subs{1},'bob')
                    obj.T= [];
                end
                varargout = cell(1,nargout);
                %varargout{1} = obj;
            else
                varargout = cell(1,nargout);
                [varargout{:}] = builtin('subsref',obj,S);
            end
             
        end
        %}
    end
    
    methods (Static)
    
        function [T] = makeT(varargin)
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
            T2 = [[cos(V(3)) sin(V(3)) 0];[-sin(V(3)) cos(V(3)) 0];[0 0 1]];
            % stretch
            T3 = zeros(3);
            T3(end,end) = 1;
            T3(1,1) = T3(1,1) + V(4);
            T3(2,2) = T3(2,2) + V(5);
            T = T1*T2*T3;
        end
        
        function [func] = makeTF(dX_range,dR_range,dS_range)
            %sgn = @(X)(X==0) + sign(X);
            glf = @(X,M,m)((abs(X)>=1)*(sign(X)+(X==0))*M) + ...
                          ((abs(X)<1)*((M-m)*X+(sign(X)+(X==0))*m));
            dX1 = @(X)glf(X,dX_range(1,1),dX_range(2,1));
            dX2 = @(X)glf(X,dX_range(1,2),dX_range(2,2));
            %dR = @(X)glf(X,dR_range(2),dR_range(1));
            dR = @(X)pi*X;
            dS1 = @(X)glf(X,dS_range(1,1),dS_range(2,1));
            dS2 = @(X)glf(X,dS_range(1,2),dS_range(2,2));
            % stack functionss
            F{1} = dX1;
            F{2} = dX2;
            F{3} = dR;
            F{4} = dS1;
            F{5} = dS2;
            % give generating function
            func = @(X)fwdT.makeT(arrayfun(@(f,X)F{f}(X),(1:numel(X)),X));
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