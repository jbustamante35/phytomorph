function [m] = p2mn(P)
    try
        %{
        %if numel(P)==5;P = P / P(end);end
        if size(P,2) == 5
            here = 1;
        end
        %}
        
        % percent less than max for width measurement
        per = .8;
        
        
        if nargin == 1
            % max value of the domain
            maxX = 2000;
            % course number of points
            xN = 1000;
            %xN = 200;
            % resolution N
            %resN = 10000000;
            % lower the resolution for speed
            resN = 10000;
            % resolution of the finest size
            dxZoom = maxX/resN;
            % zoom window size - "radius"
            widZoom = 2;
            % number of points in the zoom window to achieve the desired
            % zoom amount
            nZooom = round((2*widZoom)/dxZoom);
            % the course domain
            Xvalues = linspace(0,maxX,xN);
        end


        regrFunc = @(X,P)(P(3).*P(2).*exp(-P(3).*(X - P(1)))).*(P(4).*(exp(-P(3).*(X - P(1))) + 1).^(P(4).^-1 + 1)).^-1;


        m = zeros(size(P,1),4);
        
        for e = 1:size(P,1)
            
            tic
            % evaluate the function
            Y = regrFunc(Xvalues,P(e,:));
            % find the peak location
            [peakV,idx] = max(Y);
            % find the location in dommain
            peakX = Xvalues(idx);
            toc
            fprintf(['Done with init peak.\n']);
            
            
            tic
            p = peakX;
            F = @(X)regrFunc(X,P(e,:));
            W.numP = 10;
            W.width = widZoom;
            L.totalLoop = 0;
            L.resolutionValue = 0;
            L.loopCompression = .5;
            L.currentRes = @(J)J;
            L.halt = @(loop,res)loop < 5;
            isoNewP = @(Y)find(Y==max(Y));
            p = ifunc(p,F,W,isoNewP,L);
            toc
            fprintf(['Done with stable peak.\n']);
            
            
            
            

            
            tic
            peakZoomX = linspace(peakX-widZoom,peakX+widZoom,nZooom);
            Ypeak = regrFunc(peakZoomX,P(e,:));
            [peakV,idx] = max(Ypeak);
            peakX = peakZoomX(idx);
            toc
            fprintf(['Done with old peak.\n']);



            % find tip and base points to zoom at
            domain = Y > peakV*per;
            %tipDomain = Xvalues(find((Xvalues < peakX) & domain));
            tipDomain = Xvalues(((Xvalues < peakX) & domain));
            %baseDomain = Xvalues(find((Xvalues > peakX) & domain));
            baseDomain = Xvalues(((Xvalues > peakX) & domain));
            tipDomain = tipDomain(1);
            baseDomain = baseDomain(end);

            
            tic
            p = tipDomain;
            F = @(X)regrFunc(X,P(e,:));
            V = peakV*per;
            W.numP = 100;
            W.width = widZoom;
            L.totalLoop = 3;
            L.resolutionValue = 0;
            L.loopCompression = .5;
            L.currentRes = @(p)-log10(abs(F(p)-V));
            L.halt = @(loop,res)res<5;
            isoNewP = @(Y)find(abs((Y-V))==min(abs((Y-V))));
            stable_tipDomain = ifunc(p,F,W,isoNewP,L);
            toc
            fprintf(['Done with stable tip.\n']);
            
            tic
            % OLD ZOOM
            % zoom on tip
            tipZoomX = linspace(tipDomain-widZoom,tipDomain+widZoom,nZooom);
            Ytip = regrFunc(tipZoomX,P(e,:));
            %tipDomain = tipZoomX(find(Ytip > peakV*per));
            tipDomain = tipZoomX((Ytip > peakV*per));
            tipDomain = tipDomain(1);
            regrFunc(tipDomain,P(e,:));
            toc
            fprintf(['Done with old tip.\n']);
            
              
            tic
            p = baseDomain;
            F = @(X)regrFunc(X,P(e,:));
            V = peakV*per;
            W.numP = 100;
            W.width = widZoom;
            L.totalLoop = 3;
            L.resolutionValue = 0;
            L.loopCompression = .5;
            L.currentRes = @(p)-log10(abs(F(p)-V));
            L.halt = @(loop,res)res<5;
            isoNewP = @(Y)find(abs((Y-V))==min(abs((Y-V))));
            stable_baseDomain = ifunc(p,F,W,isoNewP,L);
            toc
            fprintf(['Done with stable base.\n']);
            
            
            
            %{
            tic;
            % zoom on tip
            baseZoomX = linspace(baseDomain-widZoom,baseDomain+widZoom,nZooom);
            Ybase = regrFunc(baseZoomX,P(e,:));
            %baseDomain = baseZoomX(find(Ybase > peakV*per));
            baseDomain = baseZoomX((Ybase > peakV*per));
            baseDomain = baseDomain(end);
            toc
            fprintf(['Done with old base.\n']);
            %}
            
            
            
            
            % make width measurement
            width = baseDomain - tipDomain;
            % make ratio measurement
            zoneR = (baseDomain-peakX)/(peakX-tipDomain);
            m(e,:) = [peakX peakV width zoneR];
        end
    catch ME
        ME;
    end
end