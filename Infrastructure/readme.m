%% Title: infrastructure
% This project will contain the functions to communicate with
% infrastructure.
%% HT Resources including squid requirements are moved to HTCresources.m
edit HTCresources
%% needs pushing me back to binary stores
% I need to know if objects are the same
%% setup dev machine to work with submit machine
% 1) make key
% 2) scp ~/.ssh/id_rsa.pub uname@submithost:.ssh/authorized_keys/
%% setup submit host for use with s3/mc/minio commands
% wget http://dl.min.io/client/mc/release/linux-amd64/mc
% chmod +x mc
% scp ~/config.json uname@submithost:.mc/config.json
% mkdir ~/bin
% mv mc ~/bin ELSE mod .bash-rc ISH thingy
%% issueTicket
% works - need to write tests
%% edit random_shell_commands.m
%% original save
% project.save('infrastructure')
project.save()
%% formal func for quinoa and example images
testF = m95('seedPlots','test');
inFile = '/mnt/spaldingdata/nate/octerineDataStore/craineBrothers/fieldQuinoa/dataSet1/IMG_9508.JPG';

testF1 = testF{inFile,inFile};
testF2 = testF1{false};

%% load infrastructure
project.load('infrastructure')
%% May 18, 2020: work on dag
%% May 19, 2020: work on resource list for jobs
FilePath = '/squid/ndmiller/';
%FilePath = '/mnt/spaldingdata/nate/test/';
[FileList] = dig(FilePath,{},{},false);
testFile = file(FileList{1});
testFile2 = file(FileList{2});
testFile3 = file(FileList{3});
%% make a resource list and add a zip file - tested loosely
rList = htcResourceList();
rList.add(testFile,'MCR');
rList.add(testFile,'MCR1');
rList.add(testFile,'MCR2');
testFileMOD = file(FileList{3});
rList.modData('MCR',testFileMOD);

%% find the matlab mcr
rFile = store.find('?,"\[required=true\]"');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% modeling a compute node as follows:
% there are energy rooms available.  each room as different energy types
% and levels.  requesting a room, bringing in resources, and setting up the
% baseplate on the energy stand in the center of the room. the base plate
% is the main.sh.  the resources are brought in via condor_file_xfer. the
% main.sh sets up the resources which may be packaged in boxes (zips).
% there maybe many requests for rooms. the submit file contains this.
%% make a collection of htcComputeNodes
l1 = htcComputeNodeCollection('Level1');
l1.generateSUBMITfile('./');
%% imread - reworked
bf = which('imread.m');
%bfNew = strrep(bf,'.p','_org.m');
bfNew = '/mnt/scratch1/phytomorph_dev/helperFunctions/imread/imread_org.m';
copyfile(bf,bfNew)
%% 
project.save()
%% one importand component is the memory profiler
% varLogger
% autoCommenter
%%
cmdText = autoCommenter.buildCMD('testFunction',4);
cmdText
cmdText_next = autoCommenter.lambdaCMD(cmdText,{'hello','world'},[2 4]);
cmdText_next

cmd = autoCommenter.buildCMD('varLogger',3);
cmd = autoCommenter.lambdaCMD(cmd,{'whos','#var1#','''f1'''});

%% to look at the tests for digs
edit digsTest.m
%% to look at the tests for digs
edit syncTest.m
%% TEST FOR NEW CFLOW
FilePath = '/squid/ndmiller/';
[squidFileList] = dig(FilePath,{},{},false);
testFile = file(squidFileList{1});

testZip = file(squidFileList{end-2});
resourceList.associatedCMD(testZip);
resourceList.add(testZip,'matlabMCR');

FilePath = '/chtc/testdata/';
[hotFileList] = dig(FilePath,{},{},false);
testFile2 = file(hotFileList{1});

resourceList = htcResourceList();
resourceList.add(testFile,'test');
curlCmds = resourceList.generateCurlCmd();

%% test for file gets - listed below are the 'scripts'
% for dig on hot or local and xfer to temp local
FileExt = {'tiff','tif','jpg','TIF','TIFF','PNG'};
FilePath = '/chtc/testb/';
FilePath = '/iplant/home/nmiller/maizeData/cobData/spaldingEars_2016/';
FileList = dig(FilePath,{},FileExt,false);
target = makeTempLocation();
xferGet(FileList,target,'clip:source',2);
%% test xfer from local to iplant
FileExt = {'tiff','tif','jpg','TIF','TIFF','PNG'};
FilePath = '/chtc/testb/';
FilePath = '/iplant/home/nmiller/JUNKY/';
FileList = dig(FilePath,{},FileExt,false);
target = makeTempLocation();
xferGet(FileList,target,'clip:source',2);
%% test xfer from local to iplant
FileExt = {'tiff','tif','jpg','TIF','TIFF','PNG'};
FilePath = '/mnt/spaldingdata/nate/octerineDataStore/testData/unitTests/';
FileList = dig(FilePath,{},FileExt,false);
target = '/iplant/home/nmiller/JUNKY/';
xferGet(FileList,target,'clip:source',2);
%% test xfer from iplant to cold
FileExt = {'tiff','tif','jpg','TIF','TIFF','PNG'};
FilePath = '/iplant/home/nmiller/JUNKY/';
FileList = dig(FilePath,{},FileExt,false);
target = '/mnt/spaldingarchive/quickTest2/';
xfer('xfer',FileList,target,'clip:source',2);
%% test the bulk transfer code(s)
tmpPath = makeTempLocation();
localList = xfer_get(FileList,tmpPath,false,true);
tmpPath = makeTempLocation();
ilocalList = xfer_get(iFileList(1:5),tmpPath,false,true);
%% test xfer file lists
target = makeTempLocation();

r = xferGet(source,target,varargin);

%%
genUpload('servicetest/yes/','10m')
%%
% test imread
I = imread(lFileList{1});
%% buld up database for ldig
% shell command to index the network drives - indexNetworkDrives.sh
%% test the autocommenter
% matlab crashed 
[newFile,newFunction] = autoCommenter.copyFunctionFile('seedPlots','_TEST');
commenter = autoCommenter();
[newFuncHandle,newFunctionName,newFunctionFile] = commenter.moniterVars(newFunction);
%%



