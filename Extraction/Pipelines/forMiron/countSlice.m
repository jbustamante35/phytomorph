function [count] = countSlice(countNet,testSlice,testMaskSlice,boxSZ)


    tmp =  testMaskSlice > .5;
    R = regionprops(tmp,'BoundingBox');
    box = R.BoundingBox;
    box(3:4) = boxSZ;
    tmpM = imcrop(testMaskSlice,box);
    testSliceToTest = [];
    for sl = 1:size(testSlice,4)
        tmp = imcrop(squeeze(testSlice(:,:,1,sl)),box);
        testSliceToTest(:,:,1,sl) = tmp.*tmpM;
    end
    u_test = squeeze(mean(testSliceToTest,4));
    szU = size(u_test);
    szU = [szU(1:2) 1 1];
    u_test = reshape(u_test,szU);
    s_test = squeeze(std(testSliceToTest,1,4));
    s_test = reshape(s_test,szU);
    z_test = zeros(size(s_test));
    x_test = cat(3,u_test,s_test,z_test);
    count = countNet.predict(x_test);
end