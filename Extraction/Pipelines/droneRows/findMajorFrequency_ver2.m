function [f,phase] = findMajorFrequency_ver2(subI,PeriodRangeFilter)
    PeriodRangeFilter = round((PeriodRangeFilter.^-1)*size(subI,2));

    % Period Range Filter - band pass
    dilateSearchRange = 31;
    filterLength = 400; % for guosheng
    filterLength = 100; % for guosheng
    % create the image to process
    toPROC = bsxfun(@minus,subI,mean(subI,2));
    % make a back ground signal
    BK = imfilter(mean(toPROC,1),fspecial('average',[1 filterLength]),'replicate');
    %toPROC = bsxfun(@minus,subI,BK);
    toPROC = bsxfun(@minus,subI,mean(subI,2));
    sig = fft(toPROC,[],2);
    %{
    WIN = 10;
    for s = 1:1:(size(toPROC,1)-WIN)
        tsig = fft(toPROC(s:(s+WIN),:),[],2);
        sig(s,:) = mean(tsig,1);
    end
    
    
    FRONT = repmat(sig(1,:),[(WIN)/2 1]);
    BACK = repmat(sig(end,:),[(WIN)/2 1]);
    sig = [FRONT;sig;BACK];
    %}
   
    sig2 = mean(abs(sig),1);
    
    
    %{
    if ~isempty(F)
        sig2 = -log(normpdf(sig2,F,FF));
    end
    %}
    
    sigRAW = sig2;
    ssig2 = imfilter(sig2,fspecial('average',[1 41]),'replicate');
    
    
    %{
    if disp
        figure;
        plot(ssig2,'k')
        hold on
        plot(sig2);
        figure;
    end
    %}
    
    % remove the trend line from the fft
    sig3 = sig2 - ssig2;
    %{
    if disp
        plot(sig3);
    end
    %}
    
    
    % search for the peak
    % change dilate amount from 51 to 31
    %cutoff = round(size(subI,2)*fCut/2);
    
    sig3(1:PeriodRangeFilter(2)) = 0;
    sig3(PeriodRangeFilter(1):end) = 0;
    
    %sig3(1:cutoff) = 0;
    sig4 = imdilate(sig3,strel('disk',dilateSearchRange,0)) == sig3;
    
    
    %{
    if isempty(fIDX)
        fidx = find(sig4);
    else
        fidx(1) = fIDX;
    end
    if disp
        hold on
        plot(sig4*30);
    end
    %}
    %fidx = find(sig4);
    [~,fidx] = max(sig3);
    f = (size(subI,2)/(fidx(1)-1))^-1;
    phase = angle(sig(:,fidx(1)));
    %{
    %fidx = fidx - 1;
    if ~isempty(fidx)
        f = (size(subI,2)/(fidx(1)-1))^-1;
        phase = angle(sig(:,fidx(1)));
    else
        f = 0;
        phase = zeros(size(subI,1),1);
    end
    %}
    
    
    
end

%{
    f = findMajorFrequency(subI,50);

%}