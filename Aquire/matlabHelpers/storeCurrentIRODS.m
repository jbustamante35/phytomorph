function [user] = storeCurrentIRODS()
    sourceFolder = '~/.irods/';
    targetFolder = '~/phytoMorphTK/.irodsStore/';
    mmkdir(targetFolder);
    json = fileread([sourceFolder 'irods_environment.json']);
    json = jsondecode(json);
    user = json.irods_user_name;
    
    sFile = [sourceFolder 'irods_environment.json'];
    tFile = [targetFolder user '_irods_environment.json'];
    copyfile(sFile,tFile);
    
    sFile = [sourceFolder 'pwfile'];
    tFile = [targetFolder user '_pwfile'];
    copyfile(sFile,tFile);
end