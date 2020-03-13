function [] = retoreIRODSuser(user)
    sourceFolder = '~/phytoMorphTK/.irodsStore/';
    targetFolder = '~/.irods/';
    
    tFile = [targetFolder 'irods_environment.json'];
    sFile = [sourceFolder user '_irods_environment.json'];
    copyfile(sFile,tFile,'f');
    
    tFile = [targetFolder 'pwfile'];
    sFile = [sourceFolder user '_pwfile'];
    
    copyfile(sFile,tFile,'f');

end