%% edit HTC resources
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test prodecure for the stucture of generating jobs
% create an object store for squid files
init = true;
if init
    % scan for files
    FilePath = '/squid/ndmiller/';
    FileList = dig(FilePath,{},{},false,true);
end
%% make the resource store
store = htcResourceStore('/mnt/scratch1/htcResources/');
%% add resources to the store
store.add(FileList{2},'[pipeline=irods]');
store.add(FileList{3},'[pipeline=compliedMatlab][mcr=840]');
store.add(FileList{4},'[pipeline=compliedMatlab][mcr=920]');
store.add(FileList{6},'[pipeline=qrcodes]');
store.add(FileList{7},'[pipeline=qrcodes]');
store.add(FileList{14},'[pipeline=dcraw]');
store.add(FileList{20},'[pipeline=lcms]');
store.add(FileList{21},'[pipeline=sqlite]');
store.add(FileList{22},'[pipeline=compliedMatlab][mcr=717]');
store.add(FileList{29},'[pipeline=compliedMatlab][mcr=930]');
%% list the resoruces
store.list()
