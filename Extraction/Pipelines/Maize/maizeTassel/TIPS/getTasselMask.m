function [tasselMask,I] = getTasselMask(I,reSZ,oPath)
    if nargin < 3;oPath = '';end
    
    % can pass image or file name
    if ischar(I);fileName = I;I = double(imread(I))/255;end
    if size(I,1) < size(I,2)
        I = imrotate(I,-90);
    end
    
    % make the binary mask
    tasselMask = thresholdTasselImage(I, reSZ);
    
    
    
    
    
    
    % write to disk if needed
    if ~isempty(oPath)
       
        
        [p,n,ex] = fileparts(fileName);
        oPath = [oPath n filesep];
        mmkdir(oPath);
        
        maskName = [oPath 'mask.tif'];
        imageName = [oPath 'image.tif'];
        overName = [oPath 'overlay.tif'];
        
        out = flattenMaskOverlay(I,logical(tasselMask));
        
        imwrite(tasselMask,maskName);
        imwrite(I,imageName);
        imwrite(out,overName);
    end
end
    