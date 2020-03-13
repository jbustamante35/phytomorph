function [] = varLogger(s,lineNum,fileID)
    %global varLog
    
    global totalMem
    global totalCount
    global lineN;
    
    %global liveChart
    
    %if isempty(liveChart);liveChart = figure;end
    
    if isempty(lineN);lineN = [];end
    
    if isempty(totalMem);totalMem = 0;end
    
    if isempty(totalCount);totalCount = 1;end
    
    if totalCount > numel(totalMem);totalMem(totalCount) = 0;end
    
    for e = 1:numel(s)
        totalMem(totalCount) = totalMem(totalCount) + s(e).bytes;
    end
    
    lineN = [lineN lineNum];
    
    
        
    %figure(liveChart);
    %plot(totalMem)
    %drawnow
    
    totalCount = totalCount + 1;
    
    %varLog(fileID) = s;
end