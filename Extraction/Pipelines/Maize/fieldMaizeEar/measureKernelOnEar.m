function [T] = measureKernelOnEar(subI,subMask,filterSize,windowSize,gridSites,CHUNK,fftFilerValue)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this function assumes that the ear is vertical
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get width profile
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    WIDTH = sum(subMask,2);
    fidx = find(WIDTH ~= 0);
    uW = mean(WIDTH(fidx));
    PS.average_WIDTH = uW;
    % store width profile
    PS.widthProfile = [interp1(1:numel(WIDTH),WIDTH,linspace(1,numel(WIDTH),1000))];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % process the grayscale image, take gradient, look at pos and neg
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if filterSize > 0
        h = fspecial('average',filterSize);            
        subI = imfilter(subI,h);
    end
    % look at gradient of grayscale image - could be yellow
    [~,subI] = gradient(subI);
    % smooth cols gradient
    subI = imfilter(single(subI),fspecial('gaussian',[41 41],11));
    % find the beginning and ending of kernels - pos and neg
    gMSK = subI > 0;
    lMSK = subI < 0;
    % strip the image border to zero
    subMask(1,:) = 0;
    subMask(:,end) = 0;
    subMask(end,:) = 0;
    % init vars to measuure image
    tG = {};
    tL = {};
    for r = 1:numel(windowSize)
        fprintf(['starting with fft window ' num2str(r) ':' num2str(numel(windowSize)) '\n']);
        % set the current window size
        dR = windowSize(r);
        % errode such that the fft window samples only ear image
        toMeasure = imerode(subMask,strel('rectangle',[2*dR+1 20]));
        % call to measure kernel period of the rising edge
        [Tg,tG{r}] = measureSingleEarFromField(gMSK.*subI,toMeasure,gridSites,dR,CHUNK,fftFilerValue);
        % call to measure the kernel period of the falling edge
        [Tl,tL{r}] = measureSingleEarFromField(lMSK.*subI,toMeasure,gridSites,dR,CHUNK,fftFilerValue);
        % stack results together
        MT(r,:) = [Tg Tl];
        % average of rising and falling edge period
        T(r) = nanmean([Tg Tl]);                
        fprintf(['ending with fft window ' num2str(r) ':' num2str(numel(windowSize)) '\n']);
    end
    fprintf(['ending with ear.\n']);
end