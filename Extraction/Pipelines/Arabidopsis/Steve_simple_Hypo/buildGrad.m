function [] = buildGrad(mag,eigVal,cIn,mE,mU,szD,featureDomains,fdSZ,disp)

    
    
    stepMAG = 1;
    
    dpath = [];
    for loop = 1:100
        dpath = [dpath cIn];
        gr = [];
        [dataOut_cp] = getDataPackage(cIn,mE,mU,szD,featureDomains,fdSZ,disp);
        dataOut_cp.features(2)
        %waitforbuttonpress
        for e = 1:numel(cIn)


            tmpC_pos = cIn;
            tmpC_pos(e) = tmpC_pos(e) + mag*eigVal(e);
            [dataOut_pos] = getDataPackage(tmpC_pos,mE,mU,szD,featureDomains,fdSZ,false);
            %tr_p = [dataOut_pos.features';dataOut_pos.s1(:);dataOut_pos.s2(:);dataOut_pos.s3(:)];
            tr_p = [dataOut_pos.features';dataOut_pos.s1(:);dataOut_pos.s3(:)];

            


            tmpC_neg = cIn;
            tmpC_neg(e) = tmpC_neg(e) - mag*eigVal(e);
            [dataOut_neg] = getDataPackage(tmpC_neg,mE,mU,szD,featureDomains,fdSZ,false);
            tr_n = [dataOut_neg.features';dataOut_neg.s1(:);dataOut_neg.s2(:);dataOut_neg.s3(:)];
            tr_n = [dataOut_neg.features';dataOut_neg.s1(:);dataOut_neg.s3(:)];


            gr = [gr ([tr_p - tr_n])/(2*mag*eigVal(e))];
            e
        end

        %[U,S,V] = svd(gr);
        %{

        dC = zeros(size(gr,2),1);
        dC(1) = 1;

        dC = V(1,:)';
        dF = (U*S*V*dC)';
        %}

        %dF = U(:,1);
        dF = zeros(size(gr,1),1);
        %dF(2) = 1;
        dF(3) = 1;
        dF(4) = 1;
        
        %dC = V*pinv(S)*U'*dF;
        dC = pinv(gr)*dF;
        %dC = dC.*eigVal;


        cIn = cIn + stepMAG*dC;
    end
    
    here = 1;
        

end