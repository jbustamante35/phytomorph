function [acc] = frameFunc(videoFile,frames,func,disp,acc)

    M = VideoReader(videoFile);

    T = M.NumberofFrames;
    HID = M.Height;   
    WID = M.Width;
    DUR = M.Duration;
    FPS = T/DUR;
    
    M = VideoReader(videoFile);

    if frames(end) == Inf;frames(end) = T;end
    if frames(1) == -Inf;frames(1) = 1;end

    % make sequence
    frames = frames(1):frames(2):frames(3);

    % convert frame sequence to seconds
    %frames = frames *(FPS)^-1;
    if nargin == 4;acc = ones(HID,WID,1);end

    
    for f = 1:numel(frames)
        tmp = M.read(frames(f));
        tmp = double(tmp)/255;
        tmp = tmp(:,:,1);
        acc = func(acc,tmp);
        if disp
            imshow(acc,[]);
            drawnow
        end
    end
end