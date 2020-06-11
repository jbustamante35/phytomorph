function [I,frames,FPS] = readSkipped(videoFile,startFrame,skipN,endFrame)

    M = VideoReader(videoFile);
    T = M.NumberofFrames;
    HID = M.Height;   
    WID = M.Width;
    DUR = M.Duration;
    FPS = T/DUR;
    
    M = VideoReader(videoFile);
    if endFrame == Inf;endFrame = T*(FPS)^-1;end
    if startFrame == -Inf;startFrame = 1;end
    frames = linspace(startFrame,endFrame,skipN);
    I = zeros(HID,WID,skipN);
    for f = 1:numel(frames)
        tmp = M.read(frames(f));
        tmp = double(tmp)/255;
        I(:,:,f) = tmp(:,:,1);
    end
end