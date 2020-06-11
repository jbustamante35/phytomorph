function [slice,szComplex] = readSlice(videoFile,F,W)

    M = VideoReader(videoFile);
    T = M.NumberofFrames;
    HID = M.Height;   
    WID = M.Width;
    DUR = M.Duration;
    FPS = T/DUR;
    szComplex = [HID,WID,T];
    M = VideoReader(videoFile);
    M.CurrentTime = (F-W)*(FPS)^-1;
    channel = 1;
    for e = 1:(2*W+1)
        tmp = M.readFrame();
        slice(:,:,channel,e) = double(tmp(:,:,1))/255;
    end
    
end