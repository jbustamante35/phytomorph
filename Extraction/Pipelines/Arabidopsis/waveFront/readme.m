%% readme
% project Name: waveFront
%%
% date: April 23, 2020
% Art and I did some work to speed up the calculation. Project is entered
% into the project manager.
% init save command: project.save('waveFront');
%% art dataset for testing
dataPath = '/mnt/tetra/nate/arthur/waveFront/';
FileExt = {'czi'};
project.attachDataCollection('Initial Dataset',dataPath,FileExt);
%%  attach dataset for single tiff stack
dataPath = '/mnt/spaldingdata/nate/octerineDataStore/waveFront/';
FileExt = {'tif','tiff'};
project.attachDataCollection('Single Test',dataPath,FileExt);
%% attach waveFront as main
func = @(in,out)waveFront(in,out);
project.attachFunction('main',func);
%% run single Test
project.runColletion('Single Test')
%% matsasugu - not yet attached
FilePath = '/home/nate/Downloads/glr3.3/glr3.3/';
FileList = {};
FileExt = {'tif','TIF'};
FileList = gdig(FilePath,FileList,FileExt,1);
%% save project
project.save();