classdef basisF < basisT
    % basis function - like a set of basis tensors
    % however mtimes is from the perspective of 
    % a wave function inner product with the observable
    methods
        
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = basisF(E,P,Z)
            if nargin == 0;P = [];Z = [];E = [];end
            if nargin == 1;P = 0;Z = [];end
            obj = obj@basisT(E,P,Z);
        end
        
        % use the function from basisT
        % this should be redone later
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [C] = mtimes(this,x)
            [C] = f(this,x);
        end
        
        % generate basis transformation to center of image
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [T] = generateCenterAffine(this)
            hsz = hsize(this);
            T = basisT.affine2D(flip(hsz,2)');
        end
        
        % half size of data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [hsz] = hsize(this)
            hsz = (size(this.E)-1)/2+1;
        end

        %{
        function [C] = prod(this,dim)
            here = 1;
        end
        %}
        
        % interpolate 2d-function
        % multi-purpose as cropping
        % x isa domain/pointset xor a box
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [res] = f(this,x,layer)
            % copy to not over write the input
            res = copy(this);
            % switch on size for crop box
            % if f*x isa crop box then crop
            if all(size(x(1).E) == [1 4])
                for e = 1:numel(res)
                    res(e).E = imcrop(res(e).E,x.E);
                end
            % else evaluate the image as function
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
    end

end