function [I,frames,FPS] = readSkippedSlices(videoFile,startFrame,skipN,endFrame,W)

    M = VideoReader(videoFile);
    T = M.NumberofFrames;
    HID = M.Height;   
    WID = M.Width;
    DUR = M.Duration;
    FPS = T/DUR;
    CH = 1;
    M = VideoReader(videoFile);
    if endFrame == Inf;endFrame = T*(FPS)^-1;end
    if startFrame == -Inf;startFrame = 1;end
    frames = linspace(startFrame,endFrame,skipN);
    I = zeros(HID,WID,skipN);
    frames(frames < W) = [];
    frames(frames > (endFrame - W)) = [];
    f = zeros(HID,WID,CH,2*W+1,numel(frames));
    
    for fr = 1:numel(frames)
        
        
        
        f(:,:,:,:,r) = readSlice(videoFile,frames(fr),W);
        
        
    end
end