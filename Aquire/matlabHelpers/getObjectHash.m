function [hash] = getObjectHash(object)
    tmp.uuid = object.uuid;
    tmp.type = object.type;
    tmp = jsonencode(tmp);
    hash = id2hash(tmp);
end




function [hash] = id2hash(id)
    cmd = ['echo -n ''' id ''' | sha256sum'];
    [~,hash] = system(cmd);
    hash = hash(1:end-4);
end