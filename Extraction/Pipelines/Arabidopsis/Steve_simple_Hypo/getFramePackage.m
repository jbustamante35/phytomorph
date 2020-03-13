function [dataOut] = getFramePackage(cIn,E,U,szD,disp)

    % sim data
    sim = PCA_BKPROJ_T(cIn,E,U);

    % extract aux1 data from sim
    mA1 = sim((end-5):(end-4),:);
    mA1 = mA1';
    
    % extract aux2 data from sim
    mA2 = sim((end-3):end,:);
    mA2 = reshape(mA2,[2 2]);
    
    % remove aux1-2 data from sim
    sim((end-5):end,:) = [];

    % reshape sim
    sim = reshape(sim,szD);
    
    
    % get in-out vecs
    vec = bsxfun(@minus,mA2,mA1);
    vec(1,:) = -vec(1,:);
    vecN = sum(vec.*vec,2).^-.5;
    vec = bsxfun(@times,vec,vecN);
    
     % get normal vecs for frame in 
    cen_nor(1) = ba_interp2(sim(:,:,end),mA1(1),mA1(2));
    cen_nor(2) = ba_interp2(sim(:,:,end-1),mA1(1),mA1(2));
    cen_nor = cen_nor / norm(cen_nor);
    cen_tanSet = [[cen_nor(2) -cen_nor(1)];-[cen_nor(2) -cen_nor(1)]];
    [~,sidx] = max(vec(1,:)*cen_tanSet');
    cen_tan = cen_tanSet(sidx,:);
    cen_nor = [cen_tan(2) -cen_tan(1)];
    cen_pt = mA1;
    cen_frame = [[cen_tan';0],[cen_nor';0],[cen_pt';1]];
    
    % get normal vecs for frame in 
    in_nor(1) = ba_interp2(sim(:,:,end),mA2(1,1),mA2(1,2));
    in_nor(2) = ba_interp2(sim(:,:,end-1),mA2(1,1),mA2(1,2));
    in_nor = in_nor / norm(in_nor);
    in_tanSet = [[in_nor(2) -in_nor(1)];-[in_nor(2) -in_nor(1)]];
    [~,sidx] = max(vec(2,:)*in_tanSet');
    in_tan = in_tanSet(sidx,:);
    in_nor = [in_tan(2) -in_tan(1)];
    in_pt = mA2(1,:);
    in_frame = [[in_tan';0],[in_nor';0],[in_pt';1]];
    
    
    
    % get normal vecs for frame in 
    out_nor(1) = ba_interp2(sim(:,:,end),mA2(2,1),mA2(2,2));
    out_nor(2) = ba_interp2(sim(:,:,end-1),mA2(2,1),mA2(2,2));
    out_nor = out_nor / norm(out_nor);
    out_tanSet = [[out_nor(2) -out_nor(1)];-[out_nor(2) -out_nor(1)]];
    [~,sidx] = max(vec(2,:)*out_tanSet');
    out_tan = out_tanSet(sidx,:);
    out_nor = [out_tan(2) -out_tan(1)];
    out_pt = mA2(2,:);
    out_frame = [[out_tan';0],[out_nor';0],[out_pt';1]];
        
        
    if disp
        imshow(sim(:,:,1),[]);
        hold on
        plot(mA2(1,1),mA2(1,2),'r.')
        plot(mA1(1),mA1(2),'g.')
        plot(mA2(2,1),mA2(2,2),'b.')
        quiver(mA2(1,1),mA2(1,2),vec(1,1),vec(1,2),vecN(1)^-1,'r')
        quiver(mA1(1),mA1(2),vec(2,1),vec(2,2),vecN(2)^-1,'b')


        plotFrame(in_frame,10);
        plotFrame(out_frame,10);
        plotFrame(cen_frame,10);

        axis([-10 50 -10 50])
    end

    dataOut.inFrame = in_frame;
    dataOut.cenFrame = cen_frame;
    dataOut.outFrame = out_frame;
    dataOut.image = sim;
    
end