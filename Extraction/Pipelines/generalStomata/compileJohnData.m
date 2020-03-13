FilePath = '/mnt/tetra/nate/STOMATA/counts_March_29_2019/';
FileList = {};
FileExt = {'csv'};
FileList = fdig(FilePath,FileList,FileExt,1);
%%
kp = [];
for e = 1:numel(FileList)
    if contains(FileList{e},'_count')
        kp(e) = 1;
    else
        kp(e) = 0;
    end
end
FileList = FileList(find(kp==1));
%%
for e = 1:numel(FileList)
    
    [pth,nm,ext] = fileparts(FileList{e});
    
    
    if mod(e,1000) == 4
        I = imread([pth filesep nm(1:end-6) '.jpg']);
        imshow(I,[]);
        drawnow
        
    end
    
    CNT(e) = csvread(FileList{e});
    OUT{e,1} = nm;
    OUT{e,2} = CNT(e);
    e
end
%%
cell2csv(['/mnt/tetra/nate/STOMATA/counts_March_29_2019/April_26_2019.csv'],OUT)
%%






















