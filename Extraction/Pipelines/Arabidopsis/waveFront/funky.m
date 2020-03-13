function [C] = funky(data,U1,E1,U2,E2)
    
dataSZ = size(data);
tmpData = reshape(data,[prod(dataSZ(1:2)) prod(dataSZ(3))]);
C = PCA_REPROJ_T(tmpData,E1,U1);
reduce1 = size(C,1);
C = reshape(C,[reduce1 dataSZ(3)]);



qSZ1 = size(C);
C = reshape(C,[prod(qSZ1(1:2)) 1]);
C = PCA_REPROJ_T(C,E2,U2);
end