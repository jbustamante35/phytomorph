function [r] = subfXfer(source,target,varargin)
    fprintf(['start: sync\n']);tic;
    fprintf(['source@' source '\n']);
    fprintf(['target@' target '\n']);
    if (isIRODS(source) && isLOCAL(target))
        cmd = ['iget ' source ' ' target];
    elseif (isLOCAL(source) && isIRODS(target))
        cmd = ['iput ' source ' ' target];
    elseif (isLOCAL(source) && isCHTC(target))
        cmd = ['mc --json cp ' source ' ' lower(target(2:end))];
    elseif (isCHTC(source) && isLOCAL(target))
        cmd = ['mc --json cp ' source(2:end) ' ' target];
    elseif (isCOLD(source) && isCHTC(target))
        % submit machine push to S3
        % have the submit machine perform the work
        cmd = ['~/bin/mc --json cp ' source ' ' lower(target(2:end))];
        cmd = ['ssh ndmiller@mir-submit.discovery.wisc.edu ''' cmd ''''];
    elseif (isCHTC(source) && isCOLD(target))
        % submit machine pulls from S3
        % have the submit machine perform the work
        cmd = ['~/bin/mc --json cp ' lower(source(2:end)) ' ' target];
        cmd = ['ssh ndmiller@mir-submit.discovery.wisc.edu ''' cmd ''''];
    elseif (isCHTC(source) && isIRODS(target))
        cmd = ['iput ' source ' ' target];
        cmd = ['ssh ndmiller@mir-submit.discovery.wisc.edu ''' cmd ''''];
    elseif (isIRODS(source) && isCOLD(target))
        cmd = ['iget ' source ' ' target];
        cmd = ['ssh ndmiller@mir-submit.discovery.wisc.edu ''' cmd ''''];
    elseif (isCOLD(source) && isIRODS(target))
        cmd = ['iput ' source ' ' target];
        cmd = ['ssh ndmiller@mir-submit.discovery.wisc.edu ''' cmd ''''];
    end
    [r,report] = system(cmd);
    report = strtrim(report);
    fprintf(['end: sync@' num2str(toc) '\n']);
end