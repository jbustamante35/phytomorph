function [P] = locatePoints(Q,curU,curSIG,zGrid,roughSpace,z_to_org)

    cur_to_zcur = @(p)bsxfun(@times,bsxfun(@minus,p,curU),curSIG.^-1);
    z_mGn = cur_to_zcur(roughSpace);
    P = [];

    parfor q = 1:size(Q,1)
        % map point to z cur
        targetQ = cur_to_zcur(Q(q,:));
        % rough pass isolation function
        disToP = @(p)sum(bsxfun(@minus,z_mGn,p).^2,2);
        % find min of func - vec list style
        isoNewP = @(d)find(d(:) == min(d(:)));
        % find the distance to the target
        disToTarget = @(X)sum(bsxfun(@minus,cur_to_zcur(p2mn(z_to_org(X),15)),targetQ).^2,2);
        % compute distance to target
        d = disToP(targetQ);
        % find the idx of the target in the rough space
        didx = isoNewP(d);
        % init the search
        %initP = gGrid(didx,:);
        initP = zGrid(didx,:);

        % number of points to use for searching
        numSearchPoints = 7;
        globalNumP = numSearchPoints*ones(1,size(zGrid,2));
        % level to find the points at
        floatLevel = 15;
        % max loop amount
        loopMax = 20;

        p = initP;
        F = disToTarget;
        P(q,:) = fminsearch(F,p);
    end
    P = z_to_org(P);
end