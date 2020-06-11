function [r] = xfer(type,source,target,varargin)
    switch type
        case 'get'
            r = xferGet(source,target,varargin{:});
        case 'put'
            r = xferGet(source,target,varargin{:});
        case 'sync'
            r = xferSync(source,target,varargin{:});
        case 'xfer'
             r = xferGet(source,target,varargin{:});
    end
end