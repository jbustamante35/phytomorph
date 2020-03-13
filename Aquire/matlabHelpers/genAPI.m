function [object] = genAPI(type,varargin)
    object = genCore(type);
    switch type
        case 'tileSeq'
            object.title = varargin{1};
            object.metaN = varargin{2};
            object.sampleN = varargin{3};
            object.usbPtr = varargin{4};
        case 'ptr'
            object.type = ['>' varargin{1}.type];
            object.refs = varargin{1};
        case 'tile'
            % set the tileType to 1)start 2)sample 3)stop 4) single
            object.type = [object.type 'Tile'];
        case 'message'
            object.fromPtr = varargin{1};
            object.toPtr = varargin{2};
            object.bodyPtr = varargin{3};
        case 'msgBody'
            flds = fields(varargin{1});
            for f = 1:numel(flds)
                object.(flds{f}) = varargin{1}.(flds{f});
            end
    end
end

function [object] = genCore(type)
    [~,uuid]=system('uuidgen');
    uuid = strtrim(uuid);
    object.genDate = datestr(now,'SSMMHHddmmYYYY');
    object.type = type;
    object.uuid = uuid;
end