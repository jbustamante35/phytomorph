function [object] = getAPIm(type,data)
    object = genCore(type);
    switch type
        case tileSeq
            flds = fields(data);
            for f = 1:numel(flds)
                object.(flds{f}) = data.(flds{f});
            end
    end
end

function [] = genCore(type)
    [~,uuid]=system('uuidgen');
    uuid = trimstr(uuid);
    object.genDate = datestr(now,'SSMMHHddmmYYYY');
    object.type = type;
    object.uuid = uuid;
end