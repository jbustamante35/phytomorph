function [] = viewImages(FilePath,per)
    if nargin == 1;per = .1;end
    M = 50;
    %% list out the images
    FileList = {};
    FileExt = {'tiff','tif','jpg','TIF','TIFF','PNG','png','JPG'};
    FileList = dig(FilePath,FileList,FileExt,1);
    %% perm
    FileList = FileList(randperm(numel(FileList)));
    
    N = min(round(per*numel(FileList)),M);
    for e = 1:N
        I = imread(FileList{e});
        imshow(I,[]);
        drawnow
        
    end
end

%{
FilePath='/mnt/spaldingdata/nate/octerineDataStore/craineBrothers/fieldQuinoa/Quinoa Seeds/2019_China-Qinghai/';
FilePath='/mnt/spaldingdata/nate/octerineDataStore/craineBrothers/fieldQuinoa/Quinoa Seeds/2019_Aus-Kununarra/';
FilePath='/mnt/spaldingdata/nate/octerineDataStore/craineBrothers/fieldQuinoa/Quinoa Seeds/2019_China-Shanxi/';

viewImages(FilePath);

%}