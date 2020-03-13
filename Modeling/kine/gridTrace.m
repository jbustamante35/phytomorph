function [] = gridTrace(initPoint,toMetricSpace,orgToZspace,zToOrgspace,dv1,dv2)
    nTrace = 30;
    
    alpha = .1;
    
    [traceMain] = subTracer(dv1,alpha,nTrace,initPoint,orgToZspace,zToOrgspace,toMetricSpace);
    
    
    traceMain_branch = arcLength(traceMain,'spec',5);
    
    
    for e = 1:size(traceMain_branch,1)
        tmpInit = traceMain_branch(e,:);
        [traceMainBranch{e}] = subTracer(dv2,alpha,nTrace,tmpInit,orgToZspace,zToOrgspace,toMetricSpace);
        fprintf(['Done tracing branch:' num2str(e) ':' num2str(size(traceMain_branch,1)) '\n'])
    end
    
    
    
end


function [traceMain] = subTracer(dv,alpha,nTrace,initPoint,orgToZspace,zToOrgspace,toMetricSpace)
    traceMainPos = orgToZspace(initPoint);
    [traceMainPos] = traceVec_ver2(toMetricSpace,traceMainPos,dv,nTrace,alpha);
    traceMainNeg = orgToZspace(initPoint);
    [traceMainNeg] = traceVec_ver2(toMetricSpace,traceMainNeg,-dv,nTrace,alpha);
    traceMainNeg(1,:) = [];
    traceMain = zToOrgspace([flip(traceMainPos,1);traceMainNeg]);
    traceMain = arcLength(traceMain,'spec',size(traceMain,1));
end