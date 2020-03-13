function [D] = gatherBaselineHistogram(imageBankLocation,N)

    FilePath = imageBankLocation;
    FileList = {};
    FileExt = {'tif'};
    FileList = fdig(FilePath,FileList,FileExt,1);

    if nargin == 1;N = numel(FileList);end
    
    N = min(N,numel(FileList))
    FileList = FileList(randperm(numel(FileList)));
    FileList = FileList(1:N);
    
    D = zeros(256,3,N);
    
    for e = 1:N
        I = imread(FileList{e});
        for k = 1:size(I,3)
            D(:,k,e) = imhist(I(:,:,k));
        end
        fprintf(['Done with image:' num2str(e) '.\n']);
    end
    D = mean(D,3);
    
end

%{

imageBankLocation ='/mnt/tetra/nate/projectData/maizeSeedling/NEF_to_tiff/';
N = 50;
D = gatherBaselineHistogram(imageBankLocation,N);

%}