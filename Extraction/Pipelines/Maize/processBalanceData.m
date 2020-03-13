FilePath = '/mnt/snapper/nate/cornellProject/finalORG/';
FileList = {};
FileExt = {'txt'};
FileList = fdig(FilePath,FileList,FileExt,1);
%%
oPath = '/mnt/snapper/nate/cornellProject/finalORG/';
oName = 'balanceTable.csv';
headerNames = {'fileName','imageType','balanceData'};
outSheet = {};
cnt = 1;
for e = 1:numel(FileList)
    if contains(FileList{e},'balanceData')
        FileList{e};
        [pth,nm,ext] = fileparts(FileList{e});
        data = readtext(FileList{e});
        data = strrep(data{1},'g/n','');
        data = strtrim(data);
        
        % process name
        nm = strrep(nm,'_balanceData','');
        fidx = strfind(nm,'_');
        type = nm((fidx(end)+1):end);
        nm = nm(1:(fidx(end)-1));
        
        outSheet{cnt,1} = nm;
        outSheet{cnt,2} = type;
        outSheet{cnt,3} = data;
        cnt = cnt + 1;
    end
end
outTable = table(outSheet(:,1),outSheet(:,2),outSheet(:,3),'VariableNames',headerNames);
writetable(outTable,[oPath oName]);
