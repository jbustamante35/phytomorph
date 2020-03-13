function [dataPoint] = extractVeinNetwork(dataPoint,oTemplate,fileName,delay)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% value for spuring the skeleton
lengthFilter = 15;
% value to measure the cost for spurs that connect branch to end
snipAmount = 31;
% amount to dilate the skeleton
skeletonDilateAmount = 3;
% amount to erode the mask
totalMaskErodeAmount = 3;
% n-hood size for threshold
nsz = 71;
% sensetivity for threshold
sen = .45;
sen = .5; % This might work better
% small hole filter
smallHoleFilter = 200;
% gap close value
% this parameter uses morphological dilatation and as a result will has
% problem with V-shapes
gapClose = 3; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all
%h1 = figure;
cnt = 1;
h1 = [];
for e = 1:numel(dataPoint)
   
    try
      

        dataPoint(e) = singleCotFromImage(dataPoint(e),sen,gapClose,lengthFilter,smallHoleFilter,snipAmount,skeletonDilateAmount,totalMaskErodeAmount,nsz);
        

%{
        % get and copy the fields
        f = fields(dataPoint(e));
        for r = 1:numel(f)
            dataPoint_out(cnt).(f{r}) = dataPoint(e).(f{r});
        end
        % set the point-set network parameters
        dataPoint_out(cnt).pointSets = pointSets;
        dataPoint_out(cnt).fileName = fileName;
        dataPoint_out(cnt).labelImage = rgbC;
        dataPoint_out(cnt).holesImage = holesImage;
        dataPoint_out(cnt).Adj = Adj;
        %}



        cnt = cnt + 1;
        if delay
            waitforbuttonpress
        end
    catch ME
        getReport(ME)
    end

end