%% attach dataset1
dataPath = '/mnt/spaldingdata/nate/octerineDataStore/craineBrothers/fieldQuinoa/dataSet1/';
project.attachDataCollection('Initial Dataset',dataPath);
%% attach dataset2
dataPath = '/mnt/spaldingdata/nate/octerineDataStore/craineBrothers/fieldQuinoa/dataSet2/';
project.attachDataCollection('Australia-2019',dataPath);
%% attach dataset3
dataPath = '/mnt/spaldingdata/nate/octerineDataStore/craineBrothers/fieldQuinoa/Quinoa Seeds/2019_China-Shanxi/';
project.attachDataCollection('China-Shanxi-2019',dataPath);
%% attach dataset4
dataPath = '/mnt/spaldingdata/nate/octerineDataStore/craineBrothers/fieldQuinoa/Quinoa Seeds/2019_China-Qinghai/';
project.attachDataCollection('China-Qinghai-2019',dataPath);
%% attach as main
func = @(in,out)seedPlots(in,out);
project.attachFunction('main',func);
%% refreseh collection
project.refreshCollection('Australia-2019');
%% number image
sz = project.collectionSize('Australia-2019');
%% test Australia collection
project.runCollection('Australia-2019');
%% test China-Shaxi collection
project.runCollection('China-Shanxi-2019');
%% test China-Qinghai collection
project.runCollection('China-Qinghai-2019');
%% test the formalFunc

%% test constructor for data collection after doid transistion
global store;
% create an obnect store
store = objectFSstore();
% generate a new collection
dataPath = '/mnt/spaldingdata/nate/octerineDataStore/craineBrothers/fieldQuinoa/dataSet2/';
newCollection = store.generate('dataCollection','Australia-2019',dataPath,{},false,true);
[fileName] = store.write(newCollection);
o = store.read(fileName);
%% profiles the function
project.profile('Australia-2019',10);
%% test run image - search by name
results = project.runDatum('Australia-2019','104');
%% test the memory profiler
% test the parts first - make sure the copy function works
% make a copy of seed plots
matFunction = which('seedPlots');
[newFile,newFunction] = autoCommenter.copyFunctionFile(matFunction,'_TEST');
fprintf(['newFile:' newFile '\n']);
fprintf(['newFunction:' newFunction '\n']);
% make mirror dim for new function
commenter = autoCommenter();
[newFunctionName,newFunctionFile,newfList] = commenter.makeMirrorDimension(newFunction);
% block check
commenter.hasVarBlock(newFile)
commenter.createComBlock(newFile)
commenter.hasVarBlock(newFile)
commenter.hasDesBlock(newFile)
%% test the function counting abilities
matFunction = which('seedPlots');
[newFile,newFunction] = autoCommenter.copyFunctionFile(matFunction,'_TEST');
[funcN,funcL] = autoCommenter.functionCount(newFile);
%% generate line map for matlab code
vec = autoCommenter.vectorize(newFile);
lineMap = autoCommenter.lineMap(vec);
commentMap = autoCommenter.commentMap(vec);
ctrlMaps = autoCommenter.ctrlMaps(vec);
line = autoCommenter.getLine(6,lineMap,vec);
UQ = unique(lineMap);
for u = 1:numel(UQ)
    autoCommenter.isLineInMap(4,lineMap,ctrlMaps(1,:));
end
%% interweave memory calls - this autocreates a mirror dimension
[newFuncHandle,newFunctionName,newFunctionFile] = commenter.moniterVars(newFunction);
%%
newFuncHandle(cur.dataCollections{1}.FileList{1},'')




