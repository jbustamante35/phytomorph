function [nslice] = resizeSlice(slice,per)
    % 
    nsliceTMP = imresize(slice(:,:,1,1,1),per);
    newSZ = [size(nsliceTMP) size(slice,3) size(slice,4) size(slice,5)];
    nslice = zeros(newSZ);
    % for each sample
    for h = 1:size(slice,5)
        fprintf(['done with sample:' num2str(h) ':' num2str(size(slice,5)) '\n']);
        % for each time-slice
        for e = 1:size(slice,4)
            % for each channel
            for c = 1:size(slice,3)
                 nslice(:,:,c,e,h) = imresize(slice(:,:,c,e,h),per);
            end
        end
    end
end