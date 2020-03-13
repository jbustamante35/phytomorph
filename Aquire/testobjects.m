global store;
% create an obnect store
store = objectFSstore();
% generate a computer
c1 = store.generate('com','computer1');
% generate usb port 1
u = store.generate('usbPort',1);
% invoke attachment
store.invoke(c1,'attachUSB',u);
% generate tile sequence
ts = store.generate('tileSequence','Title of One',u);

%%%%%%%%%%%%%%%%%%%%%%%%
%%% generate start tile
% generate start msg
strBody = store.generate('startMsg','file','folder','scanner',600,0,0);
% generate start tile
strTile = store.generate('tile');
% generate message for start tile
strMessage = store.generate('daqMessage',strTile,u,strBody);
% generate envelope
strEnvelope = store.generate('envelope',u,strMessage);
% attach envelope for start tile
store.invoke(strTile,'attachEnvelope',strEnvelope);

%%%%%%%%%%%%%%%%%%%%%%%%
%%% generate end tile
% generate start msg
endBody = store.generate('endMsg');
% generate start tile
endTile = store.generate('tile');
% generate message for start tile
endMessage = store.generate('daqMessage',endTile,u,endBody);
% generate envelope
stpEnvelope = store.generate('envelope',u,endMessage);
% attach envelope for start tile
store.invoke(endTile,'attachEnvelope',stpEnvelope);
%% test the file store with .m files
FilePath = '/mnt/scratch1/phytomorph_dev/Aquire/';
FileList = {};
FileExt = {'m'};
FileList = gdig(FilePath,FileList,FileExt,1);
for e = 1:numel(FileList)
    fileObject = store.generate('file',FileList{e});
end
%% test reading from ptr to object
testRead = store.read(dptr(stpEnvelope));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% processTemplate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FilePath = '/mnt/scratch1/phytomorph_dev/Aquire/pdfTemplates/';
FileList = {};
FileExt = {'pdf'};
FileList = gdig(FilePath,FileList,FileExt,1);
%%
height = 1;
width = @(height)height;
dx = 30;
qr_text_spacing = 10;
clear templateData;
for e = 1:numel(FileList)
    templateData(e) = processAveryTemplate(FileList{e},height,width,dx,qr_text_spacing);
end

%% make label collection
labelCol = labelCollection('CollectionTest','This is a set of test labels.');
pdfFiles = generateLabelSheets(labelCol,u,'/mnt/scratch1/phytomorph_dev/Aquire/phytoMeta/test.csv',templateData(1),'','');

%%
store.invoke(ts,'attachStartTile',strTile);
store.invoke(ts,'attachEndTile',endTile);
store.reattach(ts);
%% processTemplate

FilePath = '/mnt/scratch1/phytomorph_dev/Aquire/pdfTemplates/';
FileList = {};
FileExt = {'pdf'};
FileList = gdig(FilePath,FileList,FileExt,1);

height = .9;
width = @(height)height;
dx = 30;
qr_text_spacing = 10;

for e = 1:numel(FileList)
    templateData = processAveryTemplate(FileList{e},height,width,dx,qr_text_spacing);
end

%%



G1 = store.read(c.usbPorts(1));

T = store.read(c);
G = store.read(T.usbPorts(1));

msgBody = store.generate('startMsg','fileData','folderData','camera',400,0,0);
msg = store.generate('daqMessage','',u,msgBody);
msg.set('from',u);
%%

%%
%% from web
workerQueueConstant = parallel.pool.Constant(@parallel.pool.PollableDataQueue);
workerQueueClient = fetchOutputs(parfeval(@(x) x.Value,1,workerQueueConstant));
%%
%{
C = cell(1,4);
parfor e = 1:numel(C)
   % C{e} = e;
   C{e} = parallel.pool.PollableDataQueue();
end
%C = distributed(C);
%}
localChannel = cell(1,4);
localChannel = distributed(localChannel);
spmd
    tmp = parallel.pool.PollableDataQueue();
    localChannel(labindex) = {tmp};
    send(tmp,labindex)
end
localChannel = gather(localChannel);

localChannel = distributed(localChannel);
spmd
    getLocalPart(localChannel)
end
%%
clear future;
for e = 1:4
    future(e) = parfeval(pool,@daqThread,1);
end
%%
clear future
future = parfevalOnAll(pool,@daqThread,1);
%%
for e = 1:numel(future)
    IDX(e) = fetchOutputs(future(e))
end
%%


channel = parallel.pool.PollableDataQueue;
future = parfevalOnAll(pool,@daqThread,1,channel);
%fetchOutputs(future);
workChannels = poll(channel);


%% destroy pool
delete(gcp('nocreate'))
%% create small pool
numWorkers = 4;
if isempty(gcp('nocreate'))
    pool = parpool('local',numWorkers);
end
%% simple test (and build) here
deltaDelay = 30; % set the delay to 3 sec
publicFromWorkers = parallel.pool.PollableDataQueue;
% eval par on pool
future = parfevalOnAll(pool,@(X,N)stateCube(X,N),1,publicFromWorkers,deltaDelay);
stateCubeChannels = {};
cnt = 1;
while (publicFromWorkers.QueueLength ~= 0)
    tmpMsg = publicFromWorkers.poll();
    stateCubeChannels{cnt} = tmpMsg.body.data.channel;
    send(stateCubeChannels{cnt},clock);
    fprintf(['Send halt to stateCubes.\n']);
    cnt = cnt + 1;
end
%%
numWorkers = 2;
messageRouter(numWorkers)
%% create public channel from workers to main
publicFromWorkers = parallel.pool.PollableDataQueue;
future = parfevalOnAll(pool,@(X,N)stateCube(X,N),1,publicFromWorkers,numWorkers);
%% get IDs from workers
clear workerIDmsg
for e = 1:publicFromWorkers.QueueLength
    workerIDmsg(e) = poll(publicFromWorkers);
end
%% send msg for each worker over public channel

%%
clear privateChannelsToWorkers
for e = 1:publicFromWorkers.QueueLength
    privateChannelsToWorkers(e) = poll(publicFromWorkers);
end
%% stop workers
for e = 1:numel(privateChannelsToWorkers)
    send(privateChannelsToWorkers(e),'stop');
end




