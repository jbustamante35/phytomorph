function [m] = p2mn(P,floatLevel)
    try
        
        if nargin == 1;floatLevel = 15;end
        %{
        %if numel(P)==5;P = P / P(end);end
        if size(P,2) == 5
            here = 1;
        end
        %}
        
        % percent less than max for width measurement
        per = .8;
        
        
        %if nargin == 1
            % max value of the domain
            maxX = 2000;
            % course number of points
            xN = 3000;
            %xN = 200;
            % resolution N
            resN = 10000000;
            % lower the resolution for speed
            %resN = 10000;
            % resolution of the finest size
            dxZoom = maxX/resN;
            % zoom window size - "radius"
            widZoom = 2;
            % number of points in the zoom window to achieve the desired
            % zoom amount
            nZooom = round((2*widZoom)/dxZoom);
            % the course domain
            Xvalues = linspace(0,maxX,xN);
        %end

        method = 'new';
        
        regrFunc = @(X,P)(P(3).*P(2).*exp(-P(3).*(X - P(1)))).*(P(4).*(exp(-P(3).*(X - P(1))) + 1).^(P(4).^-1 + 1)).^-1;


        m = zeros(size(P,1),4);
        
        %floatLevel = 15;
        globalNumP = 100;
        loopMax = 40;
        %{
        % window zoom width
        widZoom = 2;
        % max value of the domain
        maxX = 2000;
        % course number of points
        xN = 3000;
        % the course domain
        Xvalues = linspace(0,maxX,xN);
        %}
        parfor e = 1:size(P,1)
            
            W = struct('numP',10,'width',widZoom);
            L = struct('totalLoop',0,'resolutionValue',0,'loopCompression',.5,...
                'updateCurrentResolution',@(J)J,'history',[]);
            
            % evaluate the function
            Y = regrFunc(Xvalues,P(e,:));
            % find the peak location
            [peakV,idx] = max(Y);
            % find the location in dommain
            peakX = Xvalues(idx);
            
            
            
            if strcmp(method,'new')
                % search for peak value and location
                p = peakX;
                F = @(X)regrFunc(X,P(e,:));
                W.numP = globalNumP;
                W.width = widZoom;
                L.totalLoop = 0;
                L.resolutionValue = 0;
                L.loopCompression = .1;
                L.updateCurrentResolution = @(p,history)-log10(abs(history - p));
                L.history = [];
                L.halt = @(loop,res)res >= floatLevel | loop > loopMax;
                isoNewP = @(Y)find(Y==max(Y));
                [peakX,retL] = ifunc(p,F,W,isoNewP,L);
                peakV = F(peakX);
            else
                peakZoomX = linspace(peakX-widZoom,peakX+widZoom,nZooom);
                Ypeak = regrFunc(peakZoomX,P(e,:));
                [peakV,idx] = max(Ypeak);
                peakX = peakZoomX(idx);
            end


            % find tip and base points to zoom at
            domain = Y > peakV*per;
            %tipDomain = Xvalues(find((Xvalues < peakX) & domain));
            tipDomain = Xvalues(((Xvalues < peakX) & domain));
            %baseDomain = Xvalues(find((Xvalues > peakX) & domain));
            baseDomain = Xvalues(((Xvalues > peakX) & domain));
            tipDomain = tipDomain(1);
            baseDomain = baseDomain(end);

            if strcmp(method,'new')
                % search for per*peakValue
                p = tipDomain;
                F = @(X)regrFunc(X,P(e,:));
                V = peakV*per;
                W.numP = globalNumP;
                W.width = widZoom;
                L.totalLoop = 0;
                L.resolutionValue = 0;
                L.loopCompression = .1;
                L.currentRes = @(p)-log10(abs(F(p)-V));
                L.halt = @(loop,res) res >= floatLevel  | loop > loopMax;
                isoNewP = @(Y)find(abs((Y-V))==min(abs((Y-V))));
                [tipDomain,retL] = ifunc(p,F,W,isoNewP,L);
            else
                tipZoomX = linspace(tipDomain-widZoom,tipDomain+widZoom,nZooom);
                Ytip = regrFunc(tipZoomX,P(e,:));
                %tipDomain = tipZoomX(find(Ytip > peakV*per));
                tipDomain = tipZoomX((Ytip > peakV*per));
                tipDomain = tipDomain(1);
            end
            
            if strcmp(method,'new')
                % search for per*peakValue
                p = baseDomain;
                F = @(X)regrFunc(X,P(e,:));
                V = peakV*per;
                W.numP = globalNumP;
                W.width = widZoom;
                L.totalLoop = 0;
                L.resolutionValue = 0;
                L.loopCompression = .1;
                L.currentRes = @(p)-log10(abs(F(p)-V));
                L.halt = @(loop,res) res >= floatLevel | loop > loopMax;
                isoNewP = @(Y)find(abs((Y-V))==min(abs((Y-V))));
                [baseDomain,retL] = ifunc(p,F,W,isoNewP,L);
            else
                
                % zoom on tip
                baseZoomX = linspace(baseDomain-widZoom,baseDomain+widZoom,nZooom);
                Ybase = regrFunc(baseZoomX,P(e,:));
                %baseDomain = baseZoomX(find(Ybase > peakV*per));
                baseDomain = baseZoomX((Ybase > peakV*per));
                baseDomain = baseDomain(end);
                
            end
            
            
            
            
            
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