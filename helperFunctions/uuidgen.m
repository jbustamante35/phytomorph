function [uuid] = uuidgen()
    [~,uuid] = system('uuidgen');
    uuid = strtrim(uuid);
end