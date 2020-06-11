function [DomainS,DomainG] = extendCurve(path,WIDTH_NUMP,PCA_RHO,WIDTH,SNIP,EXT)
    domainTranslation = [0 0];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sub-sample
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate curvilinear-domain
    %WIDTH_NUMP = 61;
    %PCA_RHO = 15;
    %WIDTH = 10;
    %SNIP = 5;
    Domain = genCurvilinearDomain(path,PCA_RHO,WIDTH,WIDTH_NUMP,[]);
    
    
    
    % extension on one side
    %EXT = 11;
    
    if EXT ~= 0

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        toFit = Domain(1:SNIP,:,:);
        dX = diff(toFit,1,1);
        dNOR = sum(dX.^2,3).^.5;
        ZdNOR = zeros(1,size(dNOR,2));
        dNOR = cat(1,ZdNOR,dNOR);
        dNOR = cumsum(dNOR,1);
        EXT_forward = -linspace(1,EXT+1,EXT);
        for d = 1:size(toFit,3)
            for w = 1:size(dNOR,2)
                b = robustfit(dNOR(:,w),toFit(:,w,d));
                mNOR(:,w,d) = b(1) + EXT_forward*b(2);
            end
        end
        mNOR = flip(mNOR,1);
        Domain = cat(1,mNOR,Domain);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        toFit = Domain((end-SNIP-1):end,:,:);
        dX = diff(toFit,1,1);
        dNOR = sum(dX.^2,3).^.5;
        ZdNOR = zeros(1,size(dNOR,2));
        dNOR = cat(1,ZdNOR,dNOR);
        dNOR = cumsum(dNOR,1);
        EXT_reverse = linspace(1,EXT+1,EXT);
        for d = 1:size(toFit,3)
            for w = 1:size(dNOR,2)
                 b = robustfit(dNOR(:,w),toFit(:,w,d));
                mNOR(:,w,d) = b(1) + (EXT_reverse+dNOR(end,w))*b(2);
            end
        end
        Domain = cat(1,Domain,mNOR);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    
    % reshape for sub-sampling
    dsz = size(Domain);
    DomainG = Domain;
    DomainG(:,:,1) = DomainG(:,:,1) + domainTranslation(1);
    DomainG(:,:,2) = DomainG(:,:,2) + domainTranslation(2);
    DomainS = reshape(Domain,[dsz(1)*dsz(2) dsz(3)]);
    DomainS = bsxfun(@plus,DomainS,domainTranslation);
end