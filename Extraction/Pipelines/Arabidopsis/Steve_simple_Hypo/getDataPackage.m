function [dataOut] = getDataPackage(cIn,E,U,szD,featureDomains,fdSZ,disp)
    box = [[0 0 1];[0 1 1];[1 1 1];[1 0 1];[0 0 0]];


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sim data
    sim = PCA_BKPROJ_T(cIn,E,U);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extract aux1 data from sim
    mA1 = sim((end-5):(end-4),:);
    mA1 = mA1';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extract aux2 data from sim
    mA2 = sim((end-3):end,:);
    mA2 = reshape(mA2,[2 2]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % remove aux1-2 data from sim
    sim((end-5):end,:) = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % reshape sim
    sim = reshape(sim,szD);

    
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get features of area,length,angle,Vin,Vout
    [frameImage] = logical(generateEdgeFrame([size(sim,1),size(sim,2)]));
    bin = sim(:,:,2) > graythresh(sim(:,:,2));
    bin_sim = sim(:,:,1) < graythresh(sim(:,:,1));




    edgebjects = frameImage & bin;
    R = regionprops(edgebjects,'Centroid');
    for e = 1:numel(R)
        mA2_bin(e,:) = R(e).Centroid;
    end
    [~,sidx] = sort(mA2_bin(:,2),'descend');
    mA2_bin = mA2_bin(sidx,:);

    



    [IO_vectors,IO_vectorsLEN] = getIOvectors(mA1,mA2);
    [IO_vectors_sim,IO_vectorsLEN_sim] = getIOvectors(mA1,mA2_bin);

    %{
    IO_vectors = IO_vectors_sim;
    mA2 = mA2_bin;
    %}

    [angle] = measureIOangle(IO_vectors(1,:),IO_vectors(2,:));
    [angle_sim] = measureIOangle(IO_vectors(1,:),IO_vectors(2,:));


    majorVector = mean(IO_vectors,1);
    majorVector = majorVector / norm(majorVector);
    majorAngle = atan2(-majorVector(2),majorVector(1));




    majorVector_sim = mean(IO_vectors_sim,1);
    majorVector_sim = majorVector_sim / norm(majorVector_sim);
    majorVector_sim = atan2(-majorVector_sim(2),majorVector_sim(1));

    
    features = [sum(bin(:)) sum(bin_sim(:)) majorAngle majorVector_sim angle angle_sim];

    %{
    features = [sum(bin(:)) angle];
    features = [sum(bin(:)) angle_sim];
    %}

    




    [cen_frame] = genFrame(sim,mA1,IO_vectors(1,:));




    [in_frame] = genFrame(sim,mA2(1,:),IO_vectors(2,:));
    [in_frame_sim] = genFrame(sim,mA2_bin(1,:),IO_vectors_sim(2,:));
    %in_frame = in_frame_sim;



    [out_frame] = genFrame(sim,mA2(2,:),IO_vectors(2,:));
    [out_frame_sim] = genFrame(sim,mA2_bin(2,:),IO_vectors_sim(2,:));
    %out_frame = out_frame_sim;
    



    
    d1 = (out_frame*featureDomains{1}')';
    d2 = (out_frame*featureDomains{2}')';
    d3 = (out_frame*featureDomains{3}')';
    s1 = squeeze(ba_interp2(sim,d1(:,1),d1(:,2)));
    s2 = squeeze(ba_interp2(sim,d2(:,1),d2(:,2)));
    s3 = squeeze(ba_interp2(sim,d3(:,1),d2(:,2)));
    s1 = reshape(s1,[fdSZ size(sim,3)]);
    s2 = reshape(s2,[fdSZ size(sim,3)]);
    s3 = reshape(s3,[fdSZ size(sim,3)]);
    


    ssXX = [size(sim,1) size(sim,2) 0];
    polyB = bsxfun(@times,box,ssXX);
    in1 = inpoly(d1(:,1:2),polyB(:,1:2));
    in2 = inpoly(d2(:,1:2),polyB(:,1:2));
    in3 = inpoly(d3(:,1:2),polyB(:,1:2));
    

        




    if disp
        
        
        imshow(sim(:,:,1),[]);
        hold on
        
     
        
        plot(mA2(1,1),mA2(1,2),'r.')
        plot(mA1(1),mA1(2),'g.')
        plot(mA2(2,1),mA2(2,2),'b.')



        quiver(mA2(1,1),mA2(1,2),IO_vectors(1,1),IO_vectors(1,2),IO_vectorsLEN(1)^-1,'r')
        quiver(mA1(1),mA1(2),IO_vectors(2,1),IO_vectors(2,2),IO_vectorsLEN(2)^-1,'b')

        quiver(mA2(1,1),mA2(1,2),IO_vectors_sim(1,1),IO_vectors_sim(1,2),IO_vectorsLEN_sim(1)^-1,'m')
        quiver(mA1(1),mA1(2),IO_vectors_sim(2,1),IO_vectors_sim(2,2),IO_vectorsLEN_sim(2)^-1,'c')

        
        plot(polyB(:,1),polyB(:,2),'r')

        plotFrame(in_frame,10);
        plotFrame(in_frame_sim,10,{'m' 'c'});
        plotFrame(out_frame,10);
        plotFrame(out_frame_sim,10,{'m' 'c'});
        plotFrame(cen_frame,10);

        axis([-10 50 -10 50])
        
        
        hold off
        drawnow
    end

    dataOut.inFrame = in_frame;
    dataOut.cenFrame = cen_frame;
    dataOut.outFrame = out_frame;
    dataOut.image = sim;
    dataOut.features = features;
    dataOut.s1 = s1;
    dataOut.s2 = s2;
    dataOut.s3 = s3;
    
    dataOut.in1 = in1;
    dataOut.in2 = in2;
    dataOut.in3 = in3;
end


function [IO_vectors,IO_vectorsLEN] = getIOvectors(centerPoint,IOpoints)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get in-out vecs
    vec = bsxfun(@minus,IOpoints,centerPoint);
    vec(1,:) = -vec(1,:);
    IO_vectorsLEN = sum(vec.*vec,2).^-.5;
    IO_vectors = bsxfun(@times,vec,IO_vectorsLEN);
end

function [angle] = measureIOangle(inVector,outVector)
    % measure angle between in and out vectors
    outFrame = [outVector;[outVector(2) -outVector(1)]];
    fMeasure = outFrame*inVector';
    angle = atan2(fMeasure(2),fMeasure(1))*180/pi;
end

function [returnFrame] = genFrame(sim,framePT,referenceVec)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get normal vecs for frame in 
    cen_nor(1) = ba_interp2(sim(:,:,end),framePT(1),framePT(2));
    cen_nor(2) = ba_interp2(sim(:,:,end-1),framePT(1),framePT(2));
    cen_nor = cen_nor / norm(cen_nor);
    cen_tanSet = [[cen_nor(2) -cen_nor(1)];-[cen_nor(2) -cen_nor(1)]];
    [~,sidx] = max(referenceVec*cen_tanSet');
    cen_tan = cen_tanSet(sidx,:);
    cen_nor = [cen_tan(2) -cen_tan(1)];
    returnFrame = [[cen_tan';0],[cen_nor';0],[framePT';1]];
end