classdef bCurve < bGraph
    
    
    properties
        length;
        strPoint;
        stpPoint;
        iFrame;
    end
    
    methods
        function [obj] = bCurve(T,X,F)
            toTransfer = false;
            if nargin >= 2
                if isa(X,'bGraph');toTransfer=true;end
            end
            obj = obj@bGraph(T);
            dX = diff(T,1,1);
            dL = sum(dX.*dX,2).^.5;
            obj.length = sum(dL);
            
            obj.strPoint = T(1,:)';
            obj.stpPoint = T(end,:)';
            
            if nargin == 2
                [iFrame,L] = bCurve.makeiFrame(T);
                obj.iFrame = bAffine(iFrame,L);
            elseif nargin == 3
                obj.iFrame = bAffine(F,1);
            end
            
            if toTransfer
                obj.namedSets = X.intersectSubset('all',obj);
            end
        end
        
        function [d] = snapToGraph(this,G)
            newStr = G.nearest(this.strPoint');
            newStp = G.nearest(this.stpPoint');
            this.E = [newStr;this.E;newStp];
            d = (norm(newStr-this.strPoint') + norm(newStr-this.strPoint'));
            this.strPoint = newStr';
            this.stpPoint = newStp';
        end
        
        function [N] = nScore(this)
            for e = 1:numel(this)
                try
                    F = this(e).iFrame.E;
                    F = F';
                    F(:,1) = F(:,1)*.5*this(e).iFrame.stretch;
                    tmpN = inv(F)*this(e).E';
                    tmpN = arcLength(tmpN','spec',50);
                    N(e) = bCurve(tmpN,[],F);
                catch ME
                    ME
                end
            end
        end
        
        function [N] = bProj(this)
            for e = 1:numel(this)
                try
                    F = this(e).iFrame.E;
                    tmpN = (F*this(e).E')';
                    N(e) = bCurve(tmpN);
                catch ME
                    ME
                end
            end
        end
        
        function [] = plot(this,cl)
            for e = 1:numel(this)
                if nargin == 1;cl = 'r';end
                plot(this(e).E(:,1),this(e).E(:,2),cl,'LineWidth',2);
                if ~isempty(this(e).iFrame)
                    this(e).iFrame.plot();
                end
            end
           
        end
        
        function [DomainG] = extendCurve(this,EXT,WIDTH,WIDTH_NUMP)
            [~,DomainG] = extendCurve(this.E(:,1:2),WIDTH_NUMP,15,WIDTH,5,EXT);
            %extendCurve(path,WIDTH_NUMP,PCA_RHO,WIDTH,SNIP,EXT)
        end
        
        function [sampleI] = sampleAlongCurve(this,EXT,WIDTH,WIDTH_NUMP,I)
            for e = 1:numel(this)
                [DomainG] = extendCurve(this(e),EXT,WIDTH,WIDTH_NUMP);
                sz = size(DomainG);
                DomainG = reshape(DomainG,[prod(sz(1:2)) sz(3)]);
                subI = ba_interp2(I,DomainG(:,2),DomainG(:,1));
                sampleI(:,:,:,e) = reshape(subI,[sz(1:2),size(subI,3)]);
            end
            %{
            imshow(I,[]);
            hold on
            plot(DomainG(:,2),DomainG(:,1),'.')
            plot(this.E(:,1),this.E(:,2),'r.');
            %}
        end
        
        function [] = arcLength(this,N)
            for e = 1:numel(this)
                tmp = arcLength(this(e).E(:,1:2),'spec',N);
                this(e).E = [tmp ones(size(tmp,1),1)];
            end
        end
        
        function [] = plotSampleGrid(this,EXT,WIDTH,WIDTH_NUMP,I)
              for e = 1:numel(this)
                [DomainG] = extendCurve(this(e),EXT,WIDTH,WIDTH_NUMP);
                sz = size(DomainG);
                DomainG = reshape(DomainG,[prod(sz(1:2)) sz(3)]);
                imshow(I,[]);
                hold on
                plot(DomainG(:,2),DomainG(:,1),'.')
                plot(this(e).E(:,1),this(e).E(:,2),'r.');
                waitforbuttonpress
              end
            
        end
        
    end
    
    
    methods (Static)
        
        function [F,L] = makeiFrame(X)
            M = .5*(X(end,:) + X(1,:));
            T = X(end,:) - X(1,:);
            N = [-T(2) T(1) T(3)];
            A = [0 0 1];
            L = norm(T);
            L = norm(T);
            N = N / norm(N);
            T = T / norm(T);
            F = [T' N' M']';
            
        end
    end
end