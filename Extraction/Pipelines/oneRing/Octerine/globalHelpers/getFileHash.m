function [ret] = fileHash(fileList,bytes)
    if ~iscell(fileList);fileList = {fileList};rFlag = 1;;end

    for e = 1:numel(fileList)
        fileName = fileList{e};
        if nargin == 2
            CMD = ['head -c' num2str(bytes) ' ''' fileName ''''];
            [r,returnHead] = system(CMD);
            CMD = ['echo -n ''' returnHead ''' | md5sum'];
            [r,returnHash] = system(CMD);
            returnHash = strrep(returnHash,'-',fileName);
        else
            CMD = ['md5sum ''' fileName ''''];
            [r,returnHash] = system(CMD);
            returnHash = trimstr(returnHash);
        end
        ret{e} = returnHash;
    end

    if rFlag;ret = ret{1};end
end