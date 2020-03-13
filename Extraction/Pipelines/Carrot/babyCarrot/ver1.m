%% scan for data
FilePath = '/mnt/spaldingdata/nate/octerineDataStore/carrotData/simonlab/plateScans/01-02-20/';
FileList = {};
FileExt = {'tif'};
tic
FileList = fdig(FilePath,FileList,FileExt,1);
toc
%% view the data
I = scanViewer(FileList);