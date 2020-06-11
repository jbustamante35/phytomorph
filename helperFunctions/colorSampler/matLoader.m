function [d] = matLoader(file,var,slice)
    try
        d = load(file,var);
        d = d.(var);
        if nargin == 3;d = d(:,slice);end
    catch ME
        getReport(ME)
    end
end