FilePath = '/mnt/snapper/nate/redCapTest/20180605Camera6_ver2/';
testFileList = {};
FileExt = {'tiff','tif','jpg','TIF','TIFF','PNG'};
testFileList = fdig(FilePath,testFileList,FileExt,1);
%%
emergenceNet = isPoppedNET;
frameNetCourse = whenPoppedNET_course;
frameNetFine = whenPoppedNET_fine;
nsz = [];
zU1 = U1;
zE1 = E1;
zL1 = L1;
zERR1 = errorL1;
zU2 = U2;
zE2 = E2;
zL2 = L2;
zERR2 = errorL2;
a = load(['./' 'cornPopperNNapp.mat']);
args = a.obj.editArgs;
kidx = kidx;
GMM = args.GMModel;
detectEmergence_ver3(testFileList,emergenceNet,frameNetCourse,frameNetFine,nsz,zU1,zE1,zL1,zERR1,zU2,zE2,zL2,zERR2,kidx,GMM)
