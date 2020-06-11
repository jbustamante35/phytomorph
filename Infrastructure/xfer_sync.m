function [r] = xfer_sync(source,target)
    fprintf(['start: sync\n']);tic;
    fprintf(['source@' source '\n']);
    fprintf(['target@' target '\n']);
    if (isIRODS(source) && isLOCAL(target))
        cmd = ['irsync -r i:' source ' ' target];
    elseif (isLOCAL(source) && isIRODS(target))
        cmd = ['irsync -r ' source ' i:' target];
    elseif (isLOCAL(source) && isCHTC(target))
        cmd = ['mc mirror ' source ' ' lower(target(2:end))];
    elseif (isCHTC(source) && isLOCAL(target))
        cmd = ['mc mirror ' lower(source(2:end)) ' ' target];
    elseif (isCOLD(source) && isCHTC(target))
        % submit machine push to S3
        % have the submit machine perform the work
        cmd = ['~/bin/mc mirror ' source ' ' lower(target(2:end))];
        cmd = ['ssh ndmiller@mir-submit.discovery.wisc.edu ''' cmd ''''];
    elseif (isCHTC(source) && isCOLD(target))
        % submit machine pulls from S3
        % have the submit machine perform the work
        cmd = ['~/bin/mc mirror ' lower(source(2:end)) ' ' target];
        cmd = ['ssh ndmiller@mir-submit.discovery.wisc.edu ''' cmd ''''];
    end
    [r,report] = system(cmd);
    report = strtrim(report);
    fprintf(['end: sync@' num2str(toc) '\n']);
end