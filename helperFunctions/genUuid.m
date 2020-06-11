function [uuid] = uuid()
    [~,uuid] = system('uuidgen');
    uuid = strtrim(uuid);
end