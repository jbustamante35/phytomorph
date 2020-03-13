function [m] = p2m(P,Xvalues)
    if nargin == 1
        Xvalues = linspace(0,2500,200);
    end


    xo = sym('xo');
    vf = sym('vf');
    k = sym('k');
    n = sym('n');
    X = sym('X');
   
    vel = sym('vel');
    regr = sym('regr');
    peak = sym('peak');
    gz = sym('gz');
   
        
    vel(X,xo,vf,k,n) = vf.*(1+exp(-k.*(X-xo))).^-(n.^-1);       % velocity
    regr(X,xo,vf,k,n) = diff(vel,X,1);                          % regr
    peak(X,xo,vf,k,n) = diff(regr,X,1);                         % where this is zero is peak
    gz(X,xo,vf,k,n) = diff(peak,X,1);                           % where this is max  == diff is zero is front

    

    parfor e = 1:size(P,1)
    
       
        %syms xo vf k n X g f
        g = sym('g');
        f = sym('f');
        h = sym('h');
        
        
        f(X) = peak(X,P(e,1),P(e,2),P(e,3),P(e,4));
        g(X) = regr(X,P(e,1),P(e,2),P(e,3),P(e,4));
        h(X) = gz(X,P(e,1),P(e,2),P(e,3),P(e,4));


        % find the peak location and value
        Y = g(Xvalues);
        [~,initMax] = max(vpa(Y));
        peakLocation = f == 0;
        peakX = vpasolve(peakLocation,X,Xvalues(initMax));
        peakV = g(peakX);
        
        
        %{
        %%%%%%%%%%% regr front and back via diff
        % find the start points for solution
        y = vpa(f(Xvalues));
        [~,maxY] = max(y);
        [~,minY] = min(y);
        regrTipLocation = vpasolve(h==0,X,Xvalues(maxY));
        regrBaseLocation = vpasolve(h==0,X,Xvalues(minY));
        %regrTipValue = g(regrTipLocation);
        %regrBaseValue = g(regrBaseLocation);
        width = regrBaseLocation - regrTipLocation;
        %}
        
        %%%%%%%%% regr front and back via percent
        per = .8;
        fidx = find(Y > per*peakV);
        regrTipLocation = vpasolve(g==per*peakV,Xvalues(fidx(1)));
        regrBaseLocation = vpasolve(g==per*peakV,X,Xvalues(fidx(end)));
        width = regrBaseLocation - regrTipLocation;
        
        % asymetric
        tipH = peakX - regrTipLocation;
        baseH = regrBaseLocation - peakX;
        zoneR = baseH / tipH;

        m(e,:) = double(vpa([peakX peakV width zoneR]));
        
        
        
    %{
    % display
    close all
    plot(Xvalues,Y,'k');hold on
    plot(peakX,peakV,'r*');
    plot(regrTipLocation,regrTipValue,'g*');
    plot(regrBaseLocation,regrBaseValue,'b*');
    plot([regrBaseLocation regrTipLocation],[peakV peakV],'m');
    %title(num2str(vpa(zoneR)));
    %}
        %if size(P,1) > 1;fprintf([num2str(e) ':' num2str(size(P,1)) '\n']);end
    end
end