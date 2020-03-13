function [xp,con] = invMap(x,X,jFunc,initP)
    options = optimset('MaxFunEvals',30000,'MaxIter',5000);
    options2 = optimoptions('fminunc','MaxFunctionEvaluations',5000,'StepTolerance',10^-20,'Display','iter','FiniteDifferenceType','central','TypicalX',ones(size(initP)),'FiniteDifferenceStepSize',[.05 .01 .05 .05]);
    for loop = 1:size(x,1)
        func = @(P)myContrast_kine1(X,P,x(loop,:));
        
        
        xp(loop,:) = fminunc(func,initP,options2);
        con(loop,1) = func(xp(loop,:));
        
        xp(loop,:) = fminsearch(func,xp(loop,:),options);
        con(loop,2) = func(xp(loop,:));
        
        xp(loop,:) = fminunc(func,xp(loop,:),options2);
        con(loop,3) = func(xp(loop,:));
        
        
        xp(loop,:) = fminsearch(func,xp(loop,:),options);
        con(loop,4) = func(xp(loop,:));
        
        con(loop,:);
        %waitforbuttonpress
        %x = particleswarm(func,4);
        %computeNewParameters(jFunc(X,xp(loop,:)),.8)
        
        %contrastF(xp(loop,:))
    end
    here = 1;
end

