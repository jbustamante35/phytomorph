%% look at example image with tips located
tI = imread('/home/nate/Downloads/18ch0041u1_imageData_kernel.tif');
imshow(tI,[]);
%% load the csv file
kernelCSVfile = '/home/nate/Downloads/16-Oct-2019_JSON_compiledData.csv';
data = readtext(kernelCSVfile);
%%
header = data(1,:);
data = data(2:end,:);
%%
fileList = unique(data(:,1));
%% page over unique filelist
for u =1:numel(fileList)
    fidx = find(strcmp(data(:,1),fileList{u}));
    subData = data(fidx,:);
    subDataU = mean(cell2mat(subData(:,9:509+499)),1);
    
    contourXU = subDataU(1:500);
    contourYU = subDataU(501:end);
    numToView = size(subData,1);
    numToView = 10;
    
    for tr = 1:numToView
        contourX = cell2mat(subData(tr,9:508));
        contourY = cell2mat(subData(tr,509:509+499));
        
        plot(contourXU,contourYU,'r')
        hold on
        plot(contourX,contourY,'b')
        hold off
        drawnow
        pause(.3)
    end
    
end
%%
for tr = 2:10
    
    contourX = cell2mat(data(tr,9:508));
    contourY = cell2mat(data(tr,509:509+499));
    plot(contourX,contourY)
    drawnow
    hold on
    waitforbuttonpress
end