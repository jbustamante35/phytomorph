function [] = spoolData(portName,spoolType,fileName,data,keyValueProps)
    global dataPool    
    port = dataPool(portName);
    switch spoolType
        case 'local'
            port.localSpool(fileName,data);
        case 'remote'
            port.localSpool(fileName,data);
        case 'buffered'  
            port.localSpool(fileName,data);
        case 'flushBuffer'
            port.flushBuffer();
    end 
end