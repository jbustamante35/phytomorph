function [] = imassDownload(source,target,fileExt,threads)
    CMD = ['mkdir -p ' target];
    system(CMD);

    if nargin == 3;threads = 4;end
    CMD = ['iquest --no-page "SELECT COLL_NAME, DATA_NAME WHERE COLL_NAME like ''' source ''''...
        ' AND DATA_NAME like ''%.' fileExt '''"'];
    CMD
    [status,text] = system(CMD);
    [text] = parseRecords(text);
    
    
    tmpF = tempname;
    fileID = fopen(tmpF,'w');
    for e = 1:numel(text)
        fileName = [text(e).COLL_NAME filesep text(e).DATA_NAME];
        fprintf(fileID,'%s\n',fileName);
    end
    
    
    CMD = ['xargs -P ' num2str(threads) ' -I % -d ''\n'' -a ' tmpF ' iget -f % ' target ];
    [status,text] = system(CMD);
    status
    text
end

%{
   
    target = '/mnt/spaldingdata/nate/tmpProject/';
    source = '/iplant/home/nmiller/cornellProject/finalORG/shipment2/kernelResults/%';
    imassDownload(source,target,'mat',12)

%}