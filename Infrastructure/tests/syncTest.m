%% test sync from local to CyVerse
sourcePath = '/mnt/spaldingdata/nate/octerineDataStore/maizeData/testData/';
targetPath = '/iplant/home/nmiller/testData/';
xfer_sync(sourcePath,targetPath);
%% test sync from CyVerse to local
sourcePath = '/iplant/home/nmiller/testData/';
targetPath = '/mnt/spaldingdata/nate/octerineDataStore/maizeData/testData/';
xfer_sync(sourcePath,targetPath);
%% test sync from local to CHTC
sourcePath = '/mnt/spaldingdata/nate/octerineDataStore/maizeData/testData/';
targetPath = '/chtc/testB/';
xfer_sync(sourcePath,targetPath);
%% test sync from local to CHTC
sourcePath = '/chtc/testB/';
targetPath = '/mnt/spaldingdata/nate/octerineDataStore/maizeData/testData1/';
xfer_sync(sourcePath,targetPath);
%% test sync from cold  to hot (CHTC)
sourcePath = '/mnt/spaldingarchive/jpatton_testing/';
targetPath = '/chtc/testb/';
xfer_sync(sourcePath,targetPath);
%% test sync from hot to cold (CHTC)
sourcePath = '/chtc/testb/';
targetPath = '/mnt/spaldingarchive/quickTest/';
xfer_sync(sourcePath,targetPath);