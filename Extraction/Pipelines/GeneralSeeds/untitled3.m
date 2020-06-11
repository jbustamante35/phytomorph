FilePath = '/mnt/spaldingdata/nate/octerineDataStore/craineBrothers/ColorAnalysis/';
FileSet = {};
FileExt = {'tiff','tif','jpg','TIF','TIFF','PNG'};
FileSet = sfdig(FilePath,FileSet,FileExt,1);
%% view the data
close all
for s = 1:numel(FileSet)
    FileList = FileSet{s};
    for e = 1:numel(FileList)
        I = imread(FileList{e});
        imshow(I,[]);
        drawnow
    end
end