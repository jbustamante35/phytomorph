function [kvp] = pargin(varargin)
    if mod(numel(varargin),2) ~= 0;kvp = '';return;end
    for e = 1:(numel(varargin)/2)
        idx = (e-1)*2 + 1;
        kvp(e).key = varargin{idx};
        kvp(e).value = varargin{idx+1};
    end
end

%{
    kvp = pargin('hello,'world');
%}