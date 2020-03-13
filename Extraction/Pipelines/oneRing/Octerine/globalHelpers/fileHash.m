function [ret] = fileHash(fileList,tableOutput,bytes)
    rFlag = 0;
    if nargin == 1;tableOutput = false;end
    if ~iscell(fileList);fileList = {fileList};rFlag = 1;end

    for e = 1:numel(fileList)
        fileName = fileList{e};
        if nargin == 3
            CMD = ['head -c' num2str(bytes) ' ''' fileName ''''];
            [r,returnHead] = system(CMD);
            CMD = ['echo -n ''' returnHead ''' | md5sum'];
            [r,returnHash] = system(CMD);
            returnHash = strrep(returnHash,'-',fileName);
        else
            CMD = ['md5sum ''' fileName ''''];
            [r,returnHash] = system(CMD);
            sidx = strfind(returnHash,'  ');
            returnHash(sidx(1)) = '*';
            returnHash(sidx(1)+1) = [];
            returnHash = strtrim(returnHash);
        end

        ret{e} = returnHash;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if tableOutput
        for e = 1:numel(ret)
            tmp = ret{e};
            sidx = strfind(tmp,'*');
            hashVar{e} = tmp(1:(sidx(1)-1));
            wholeFile = tmp((sidx(1)+1):end);
            [~,nm] = fileparts(wholeFile);
            wholeVar{e} = wholeFile;
            nameVar{e} = nm;
        end
        ret = table(wholeVar',nameVar',hashVar','VariableNames',{'wholeName','fileName','hashValue'});
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if rFlag;ret = ret{1};end
end