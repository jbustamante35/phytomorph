function [fileName] = getObjectFile(object)
    objectHome = getenv('PHYTO_OBJECTS');
    hash = getObjectHash(object);
    type = object.type;
    fileName = [objectHome type filesep hash(1:2) filesep hash];
end