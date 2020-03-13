function [nC] = realignCenters_singulated(fileName,C,BOX)
    for e = 1:size(C,1)
        I = getKernelImage(fileName,C(e,:),BOX);
        
         % make gray scale color image
        if size(I,3) == 1
            I = cat(3,I,I,I);
        end
        
        
        nC(e,:) = getNewCenter_singulated(I,10^3);
        nC(e,:) = (fliplr(nC(e,:))) + C(e,:);
        fprintf(['Done realigning center ' num2str(e) ':' num2str(size(C,1)) '\n']);
    end
end