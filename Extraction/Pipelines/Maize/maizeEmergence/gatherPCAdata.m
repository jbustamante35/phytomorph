function [c,e,simc] = gatherPCAdata(fileName,imageScale,sampleFunc,U,E)
    tmp = imread(fileName);
    tmp = imresize(tmp,imageScale);
    [colorSample] = sampleFunc(tmp);
    c = PCA_REPROJ_T(colorSample,E,U);
    simc = PCA_BKPROJ_T(c,E,U);
    e = sum((colorSample - simc).*(colorSample - simc),1);
end