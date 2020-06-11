%% look for image files
FilePath = '/iplant/home/nmiller/maizeData/cobData/spaldingEars_2016/';
FileList = {};
FileExt = {'tiff','tif','jpg','TIF','TIFF','PNG'};
FileList = idig(FilePath,FileList,FileExt);
%% spin out and up DAG
% setup parameters for algorithm
noe = 3;
oPath = './output/';
rPath = [];
rawImage_scaleFactor = 1;
checkBlue_scaleFactor = 1;
defaultAreaPix = 10^6;
rho = 300;
addcut = 150;
baselineBlue = 600;
colRange1 = 70;
colRange2 = 166;  
fill = 50;
toSave = 1;
toDisplay = 1;
%% try out a few local
for e = 1:numel(FileList)
    fileName = FileList{e};
    singleCobImage(fileName,noe,oPath,rPath,rawImage_scaleFactor,checkBlue_scaleFactor,defaultAreaPix,rho,addcut,baselineBlue,colRange1,colRange2,fill,toSave,toDisplay)
end 
%%
FileList = issueBulkTicket(FileList);
%% try out some on condor
func = cFlow('singleCobImage');
func.setMCRversion('v930');
for e = 1:5%numel(FileList)
    fileName = FileList{e};
    res{e} = func(fileName,noe,oPath,rPath,rawImage_scaleFactor,checkBlue_scaleFactor,defaultAreaPix,rho,addcut,baselineBlue,colRange1,colRange2,fill,toSave,toDisplay)
end
auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
func.submitDag(auth,150,150);
%%








