function [I] = scanViewer(FileList,N,skip)
    if nargin == 1;N = numel(FileList);skip = 1;end
    for e = 1:skip:N
        I = imread(FileList{e});
        R = imfinfo(FileList{e});
        R.BitDepth
        imshow(I,[]);
        if R.BitDepth ~= 1
            waitforbuttonpress
        end
        drawnow
    end
end