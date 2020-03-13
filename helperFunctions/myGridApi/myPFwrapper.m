function [varargout] = myPFwrapper(functionName,varargin)

    options = weboptions('Timeout',60);
    webPath = './';
    fprintf(['start:pulling ' functionName '.mat\n']);
    websave([webPath functionName '.mat'],['https://data.cyverse.org/dav-anon/iplant/home/nmiller/publicData/' functionName '.mat'],options);
    load([webPath functionName '.mat']);
    funcToCall = obj.func;
    delete([webPath functionName '.mat']);
    funcToCall(varargin{:});
end