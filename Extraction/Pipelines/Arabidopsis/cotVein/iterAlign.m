function [mdata,idata,Trans] = iterAlign(mdata,mtemplate,idata,RD,L)
    dispO = false;
    
    [optimizer, metric] = imregconfig('multimodal');
    
    optimizer.InitialRadius = optimizer.InitialRadius*.05;
    optimizer.Epsilon = 1.5e-4;
    optimizer.GrowthFactor = 1 + (optimizer.GrowthFactor-1)*.5;
    optimizer.MaximumIterations = 300;
    
    for loop = 1:L
        Trans{loop} = imregtform(mdata,mtemplate,'similarity',optimizer,metric,'DisplayOptimization',dispO,'PyramidLevels',4);
        mdata = imwarp(mdata,Trans{loop},'OutputView',RD);
        idata = imwarp(idata,Trans{loop},'OutputView',RD);
    end
end