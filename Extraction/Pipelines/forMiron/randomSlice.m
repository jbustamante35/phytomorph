function [f,s,F] = randomSlice(videoFiles,W,N)
    % for N samples
    for r = 1:N
        fprintf(['Sampling number:' num2str(r) ':' num2str(N) '\n']);
        % generate random frame to read
        s(r) = randi(numel(videoFiles),1);
        % select "this" video file
        videoFile = videoFiles{s(r)};
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        M = VideoReader(videoFile);
        T = M.NumberofFrames;
        HID = M.Height;   
        WID = M.Width;
        DUR = M.Duration;
        FPS = T/DUR;
        CH = 1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if r == 1;f = zeros(HID,WID,CH,2*W+1,N);end
        upper = T - (W+1);
        lower = (W+1);
        % pre-allocate the slice
        fr = round((upper - lower)*rand(1) + lower);
        % read the slice
        f(:,:,:,:,r) = readSlice(videoFile,fr,W);
        % return the frame choosen
        F(r) = fr;
    end
end