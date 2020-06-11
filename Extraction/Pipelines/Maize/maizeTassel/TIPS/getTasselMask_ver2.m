function [tasselMask,I] = getTasselMask_ver2(fileName,reSZ,AImethods,oPath)
    if nargin < 4;oPath = '';end
    [p,nm,ext] = fileparts(fileName);
    % run original method
    [tasselMask,I] = getTasselMask(fileName,reSZ);
    
    % stack data for AImethods
    szI = size(I);
    I = reshape(I,[prod(szI(1:2)) szI(3)]);
    
    
    
    for m = 1:numel(AImethods)
        fprintf(['start: ai-method:' num2str(m) ':' num2str(numel(AImethods)) '\n']);tic
        y = AImethods(m).operate(I*255);
        if size(y,2) > size(y,1);y = y';end
        if size(y,2) > 1;y = y(:,2) > .5;end
        y = reshape(y,szI(1:2));
        % not always the largest
        y = bwlarge(y);
        
        %y = imclearborder(y);
        %y = bwlarge(y);
        
        tasselMask = cat(3,tasselMask,y);
        fprintf(['start: ai-method:' num2str(m) ':' num2str(numel(AImethods)) ':' num2str(toc) '\n']);
    end
    
    if ~isempty(oPath)
        for m = 1:size(tasselMask,3)
            fileName = [oPath '{originalName_' nm '}{maskNumber_' num2str(m) '}.tif'];
            imwrite(tasselMask(:,:,m),fileName);
        end
    end
    
    % restore I and return
    I = reshape(I,szI);
    
end
    