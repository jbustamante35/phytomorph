function [block] = generateImageBlocks(blockSZ,N,state,colorList)
    block = [];
    for e = 1:N
        tmpBlock = makeColorBlock(blockSZ,colorList(state(e),:));
        block = [block,tmpBlock];
    end
end

%{

    blockSize = [100 100];
    N = 4;
    state = [1 1 1 1];
    colorList = [[1 0 0];[0 1 0];[0 0 1]];
    test = generateImageBlocks(blockSize,N,state,colorList);
    imshow(test,[]);

    close all
    for e = 1:10
        state = randi([1 3],1,N);
        test = generateImageBlocks(blockSize,N,state,colorList);
        imshow(test,[]);
        drawnow
    end


%}


