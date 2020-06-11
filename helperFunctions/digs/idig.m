function [FileList] = idig(FilePath,FileList,FileExt,verbose)    
    cmd = ['iquest --no-page "SELECT DATA_ID, DATA_NAME, COLL_NAME WHERE COLL_NAME like ''' FilePath '%''"'];
    
    
    [status,data] = system(cmd);
    if ~strcmp(FilePath,'/')
        FilePath = [FilePath '/'];
    end
    data = parseRecords(data);
    
    for e = 1:numel(data)
        [~,~,ext] = fileparts(data(e).DATA_NAME);
        ext(1) = [];
        if contains(ext,FileExt)
            FileList{end+1} = [data(e).COLL_NAME '/' data(e).DATA_NAME];
        end
    end
    
end

%{

    FilePath = '/iplant/home/nmiller/maizeData/cobData/spaldingEars_2016/';
    FileList = {};
    FileExt = {'tiff','tif','jpg','TIF','TIFF','PNG'};
    FileList = idig(FilePath,FileList,FileExt);
   

    

%}