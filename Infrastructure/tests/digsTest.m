%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dig and construct file for LOCAL
FileExt = {'tiff','tif','jpg','TIF','TIFF','PNG','zip'};
FilePath = '/mnt/spaldingdata/nate/octerineDataStore/testData/unitTests/';
[FileList] = dig(FilePath,{},FileExt,false);
testFile = file(FileList{1});
%% dig and construct file for IRODS
FilePath = '/iplant/home/nmiller/maizeData/cobData/spaldingEars_2016/';
[FileList] = dig(FilePath,{},FileExt,false);
testFile = file(FileList{1});
%% dig and construct file for HOT
FilePath = '/chtc/';
[FileList] = dig(FilePath,{},FileExt,false);
testFile = file(FileList{1});
%% dig and construct file for COLD/CHTC
FilePath = '/mnt/spaldingarchive/quickTest/';
[FileList] = dig(FilePath,{},FileExt,false);
testFile = file(FileList{1});
%% dig and construct file for SQUID
FilePath = '/squid/ndmiller/';
[FileList] = dig(FilePath,{},FileExt,false);
testFile = file(FileList{1});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% list of digs
% gdig - search using ls/dir
% fdig - search using find
% ldig - dig for files using locate
% idig - search for files at CyVerse
% wdig - searcg CHTC
%% test digs
[location] = getLocation(locationString);
%% test gdig
FilePath = '/mnt/spaldingdata/nate/octerineDataStore/testData/unitTests/';
lFileList = {};
FileExt = {'JPG'};
lFileList = gdig(FilePath,lFileList,FileExt,true);
location = getLocation(FilePath)
%% test fdig
FilePath = '/mnt/spaldingdata/nate/octerineDataStore/testData/unitTests/';
FileList = {};
FileExt = {'JPG'};
FileList = fdig(FilePath,FileList,FileExt,true);
location = getLocation(FilePath)
%% test ldig
FilePath = '/mnt/spaldingdata/nate/octerineDataStore/testData/unitTests/';
FileList = {};
FileExt = {'JPG','jpg'};
[FileList] = ldig(FilePath,FileList,FileExt,true);
location = getLocation(FilePath)
%% test idig
FilePath = '/iplant/home/nmiller/maizeData/cobData/spaldingEars_2016/';
iFileList = {};
FileExt = {'tiff','tif','jpg','TIF','TIFF','PNG'};
iFileList = idig(FilePath,iFileList,FileExt);
location = getLocation(FilePath)
%% test wdig
FilePath = '/chtc/';
FileList = {};
FileExt = {'tiff','tif','jpg','TIF','TIFF','PNG'};
FileList = wdig(FilePath,FileList,FileExt);
location = getLocation(FilePath)
%% test cold dig using find
FilePath = '/mnt/spaldingarchive/jpatton_testing/';
FileList = {};
FileExt = {'dat'};
FileList = fdig(FilePath,FileList,FileExt,true);
location = getLocation(FilePath)