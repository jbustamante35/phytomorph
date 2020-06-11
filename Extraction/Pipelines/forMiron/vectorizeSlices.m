function [v] = vectorizeSlices(slice,szV,skV,func)
    v = [];
    % loop over each sample
    for sam = 1:size(slice,5)
        fprintf(['Vectorizing sample:' num2str(sam) ':' num2str(size(slice,5)) '\n']);
        for sl = 1:size(slice,4)
            % look over each slice
            for c = 1:size(slice,3)
                tmp = im2colF(slice(:,:,c,sl,sam),szV,skV);
                % apply the func
                tmpf = func(tmp);
                if isempty(v);v = zeros([size(tmpf) size(slice,3) size(slice,4) size(slice,5)]);end 
                v(:,:,c,sl,sam) = tmpf;
            end
        end
    end
end