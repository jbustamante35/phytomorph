classdef funcP < fwdT
    
    properties
        F;
        noCache = true;
    end
    
    methods
        function [obj] = funcP(P,fT)
            obj = obj@fwdT(P,fT);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sample the function
        % x - domain
        % sz - size of domain
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [F] = getF(obj,globData,x,sz,layer,id)
            if nargin == 4;layer = 1:size(globData,3);end
            if nargin < 6;id = 1;end
            F = [];
            for e = 1:numel(obj)
                if isempty(obj(e).F) || (obj.noCache)
                    if ~isempty(sz);patchSZ = [sz numel(layer) numel(id)];end
                    T = obj(e).getT(globData);
                    x = mtimesx(T,x,'T')';
                    v = zeros([size(x,1) numel(layer) numel(id)]);
                    for i = 1:numel(id)
                        
                        v(:,:,i) = squeeze(ba_interp2(globData(:,:,:,id(i)),x(:,1),x(:,2)));
                        
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
       
        function [d] = dY(A,B,globDataA,globDataB,x,sz)
                a = A.getF(globDataA,x,sz);
                b = B.getF(globDataB,x,sz);
                d = norm(a(:)-b(:));
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