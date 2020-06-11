function [] = mmkdir(systemPath)
    if isLOCAL(systemPath)
        CMD = ['mkdir -p ''' systemPath ''''];
    elseif isIRODS(systemPath)
        CMD = ['imkdir -p ''' systemPath ''''];
    elseif isCOLD(systemPath)
        CMD = ['mkdir -p ''' systemPath ''''];
        CMD = ['ssh ndmiller@mir-submit.discovery.wisc.edu ''' CMD ''''];
    end
    system(CMD);
end