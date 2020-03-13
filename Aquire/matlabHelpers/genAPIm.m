function [object] = getAPI(type,varargin)
    object = genCore(type);
    switch type
        case 'tileSeq'
            object.title = varargin{1};
            object.metaN = varargin{2};
            object.sampleN = varargin{3};
            object.usbPtr = varargin{4};
    end
end

function [object] = genCore(type)
    [~,uuid]=system('uuidgen');
    uuid = strtrim(uuid);
    object.genDate = datestr(now,'SSMMHHddmmYYYY');
    object.type = type;
    object.uuid = uuid;
end