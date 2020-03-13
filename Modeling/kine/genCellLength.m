function [initP] = genCellLength(velFunc,regrFunc,X,para,initP,dt,TAU,cellProduction)
    % initP - init position
    % para - parameters for vel
    % dt = time step
    % TAU - total time
    % dt - in hours
    % TAU - in hours
    func = @(X,xo,vf,k,n)vf.*(1+exp(-k.*(X-xo))).^-(n.^-1);
    
    nT = round(TAU / dt);
    tm = linspace(0,TAU,nT);
    
    

    conFPS = 30^-1;
    conSPH = 60*60;
    conPPM = 1463^-1;
    
    conFPS = 1;
    conSPH = 1;
    conPPM = 1;
    
    
    
    %vel = velFunc(X,para);
    %vel = func(X,para(1),para(2),para(3),para(4));
    
    
    
    
    %vel = vel*conFPS*conSPH*conPPM;
    
    %regr = gradient(vel,conPPM);
    %regr = gradient(vel,mean(diff(X)));
    
   
    
    X = X*conPPM;
    
    
    
    for t = 1:numel(tm)
        %for p = 1:size(initP,2)
        
           
        
            p = initP(end,1);
            dp = initP(end,2);
            
            %velAtP = interp1(X,vel,initP(end,1));
            %regrAtP = interp1(X,regr,initP(end,1));
            
            velAtP = velFunc(p,para);
            regrAtP = regrFunc(p,para);
            frac = interp1(X,cellProduction,p,'spline');
            %frac = frac * dp;
            
            tmpP(1) = p + velAtP*dt;
            tmpP(2) = frac*(1+(regrAtP*dt))*dp;
            
            
            
            initP = [initP;tmpP];
            
            %{
            plot(initP(:,2))
            waitforbuttonpress
            %}
    end
    tm = [0;tm'];
    initP = [initP tm];
   %ß L = initP(:,2) - initP(:,1);
end