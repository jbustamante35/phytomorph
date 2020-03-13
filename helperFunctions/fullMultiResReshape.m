function [rI] = fullMultiResReshape(I,res,BOX)
    PAD = max((BOX*res')/2,[],2))
    for e = 1:numel(res)
        %J = imresize(I,res(e));
        tmpBOX = round(BOX*res(e));
        sI = im2colF(double(I),tmpBOX,[1 1]);
        sI = reshape(sI,[tmpBOX 1 size(sI,2)]);
        for n = 1:size(sI,4)
            rI(:,:,:,n,e) = imresize(sI(:,:,:,n),BOX);
            size(sI,4)
        end
    end
       
end