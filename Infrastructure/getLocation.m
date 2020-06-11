function [location] = getLocation(locationString)
    if isLOCAL(locationString);location = 'local';end
    if isCHTC(locationString);location = 'chtc';end
    if isIRODS(locationString);location = 'irods';end
    if isCOLD(locationString);location = 'cold';end
    if isSQUID(locationString);location = 'squid';end
end