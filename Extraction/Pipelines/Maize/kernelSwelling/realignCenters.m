function [nC] = realignCenters(fileName,C,BOX)
    debug = false;
    for e = 1:size(C,1)
        I = getKernelImage(fileName,C(e,:),BOX);
        if debug
            imshow(I,[]);
            waitforbuttonpress
        end
        nC(e,:) = getNewCenter(I,10^3);
        nC(e,:) = fliplr(nC(e,:)) + C(e,:);
        fprintf(['Done realigning center ' num2str(e) ':' num2str(size(C,1)) '\n']);
    end
end