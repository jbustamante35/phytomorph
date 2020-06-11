function [I] = randomVIslice(videoFile,W,N)
    videoFiles{1} = videoFile;
    [f,s,F] = randomSlice(videoFiles,W,1);
    %opticFlow1 = opticalFlowLK();
    %opticFlow1 = opticalFlowLKDoG();
    opticFlow1 = opticalFlowFarneback('NeighborhoodSize',N);
    for e = 1:W
        flowField = estimateFlow(opticFlow1,f(:,:,e));
    end
    %I = cat(3,f(:,:,end,1),flowField.Vx,flowField.Vy);
    I = cat(3,f(:,:,end,1),flowField.Magnitude);
end