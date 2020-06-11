function [U,E,L] = PCA_FIT_FULL_T2ws(X1,X2,COM,colorFunc)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % performs PCA (principal component analysis) using eigenvector
    % decomposition, backprojects to simulate the data and calculates the
    % error
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUTS:   
    %           M       : = data matrix
    %           COM     : = number of vectors 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUTS:  
    %           S       : = simulated signal (reconstruction)
    %           C       : = components ("unique" fingerprint)
    %           U       : = mean of data
    %           E       : = basis vectors 
    %           L       : = eig values
    %           ERR     : = error in reconstruction
    %           LAM     : = percent explained
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sigma default
    if nargin == 3
        sigma = 'largestabs';
    end
    %{
    a = functions(colorFunc);
    maskIDX = a.workspace{1}.w.workspace{1}.kidx;
    mask = zeros([226 226]);
    mask(maskIDX) = 1;
    mask = imresize(mask,[301 301]);
    imshow(mask,[]);
    %}
    
    initCOM = 20;
    [U1,E1,L1] = PCA_FIT_FULL_Tws(X1,initCOM);
    [U2,E2,L2] = PCA_FIT_FULL_Tws(X2,initCOM);
    S1 = PCA_REPROJ_T(X1,E1,U1);
    S2 = PCA_REPROJ_T(X2,E2,U2);
    
    
    for e = 1:size(X2)
        tmpI = colorFunc(X2(:,e));
        imshow(bindVec(tmpI),[]);
        title(num2str(e))
        %waitforbuttonpress
        drawnow
    end
    % view eigen vectors
    for e = 1:size(E2,2)
        tmpI = colorFunc(E2(:,e));
        imshow(bindVec(tmpI),[]);
        waitforbuttonpress
    end
   
    [sweepD] = sweepPCA(S1',E1,U1',diag(L1).^.5,1:3,5);
    [sweepD] = sweepPCA(S2',E2,U2',diag(L2).^.5,1:3,5);
    for com = 1:size(sweepD,1)
        for val = 1:size(sweepD,2)
            sig = squeeze(sweepD(com,val,:));
            tmpI = colorFunc(sig);
            imshow(tmpI/255,[0 1]);
            title([num2str(com) '-' num2str(val)]);
            drawnow
            %waitforbuttonpress
            pause(.5)
        end
    end
    
    
    toOpDim = 2;
    toKeep = min(size(S1,toOpDim),size(S2,toOpDim));
    S1 = S1(:,1:toKeep);
    S2 = S2(:,1:toKeep);
    
    [cs1,cs2,r,q1,q2,stats] = canoncorr(S1',S2');
    
    [U1,E1,L1] = PCA_FIT_FULL_Tws(cs1(:,1:3),initCOM);
    [U2,E2,L2] = PCA_FIT_FULL_Tws(cs2(:,1:3),initCOM);
    
    % take the mean
    fprintf(['PCA:start:taking off mean \n']);
    toOpDim = 2;
    U1 = mean(X1,toOpDim);
    U2 = mean(X2,toOpDim);
    X1 = bsxfun(@minus,X1,U1);
    X2 = bsxfun(@minus,X2,U2);
    fprintf(['PCA:end:taking off mean \n']);
    
    
    
    % look at covariance
    fprintf(['PCA:start:creating COV \n']);
    try
        fprintf(['i am speed.l.mcqueen \n']);
        COV12 = mtimesx(X1,X2,'T','speed');
        fprintf(['and you know that.l.mcqueen \n']);
    catch
        COV12 = X1*X2';
    end
    COV12 = COV12 / size(X1,2);
    fprintf(['PCA:end:creating COV \n']);
    % eig vector decomp
    fprintf(['PCA:start:decomposing COV \n']);
    [E12,L12] = eigs(COV12,COM(2));
   
    L12 = diag(L12);
    pidx = find(imfill(~(imag(L12) == 0),1,8) == (imag(L12) == 0));
    E12 = E12(:,pidx);
    COM(2) = numel(pidx);
    fprintf(['PCA:end:decomposing COV \n']);
    L12 = diag(L12);
    L12 = L12(pidx);
    % project X1 into common
    sX1 = mtimesx(E12,'T',X1);
    simX1 = mtimesx(E12,sX1);
    norX1 = X1 - simX1;
    % project X1 into common
    sX2 = mtimesx(E12,'T',X2);
    simX2 = mtimesx(E12,sX2);
    norX2 = X2 - simX2;
    [subU1,subE1,subL1] = PCA_FIT_FULLws(norX1,COM(1));
    [subU2,subE2,subL2] = PCA_FIT_FULLws(norX2,COM(3));
end








