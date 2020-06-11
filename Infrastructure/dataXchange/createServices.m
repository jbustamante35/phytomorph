%% make cyverse,cold store, and local store folder structure
% structureis:
% /phytoMorphService/service/in/#collaboratorID#/#siteID#/
% /phytoMorphService/collaborator/
% (user,service) pair for each location
base{1} = '/iplant/home/nmiller/';
base{2} = '/mnt/spaldingdata/';
base{3} = '/mnt/spaldingarchive/';
userLocation = 'phytoMorphService/collaborator/';
serviceLocation ='phytoMorphService/service/in/';
serviceLocation ='phytoMorphService/service/out/';
for e = 1:numel(base)
    folderName = [base{e} userLocation];
    mmkdir(folderName)
end
for e = 1:numel(base)
    folderName = [base{e} serviceLocation];
    mmkdir(folderName)
end
%% create the bucket "phytomorphservice"
bucketName = 'phytomorphservice/'
cmd = ['mc mb chtc/' bucketName];
[o.result] = system(cmd);