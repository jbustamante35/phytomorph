%% attach test stacks
dataPath = '/mnt/spaldingdata/nate/octerineDataStore/hypocotylSegments/';
project.attachDataCollection('testStacks',dataPath,'',true);
%% attach vectorized list of images
dataPath = '/mnt/spaldingdata/nate/octerineDataStore/hypocotylSegments/';
project.attachDataCollection('vectorizedTestStacks',dataPath,'',false);
%% remove collections
project.removeCollection('vectorizedTestStacks');
project.removeCollection('testStacks');
%% refresh collections
project.refreshCollection('all');
%% attach as main
func = @(in,out)op0_liteDE(in,out,false);
project.attachFunction('main',func);
%% attach single image sub function
func = @(in,out)opOnSingleImage(in,false);
project.attachFunction('processSingleImage',func);
%%
project.profile('testStacks',3);
%%
project.profile('vectorizedTestStacks',10,'processSingleImage');
%%
[L,WID,dB,P,loc] = opOnSingleImage(imageName,disp);
%% profile on test stacks
project.profile('testStacks',3);
%% 