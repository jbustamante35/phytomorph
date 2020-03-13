function [str] = getValidRead(prompt)
    str = '';
    while ~validPSON(str)
        str = input(prompt,'s');
        % if valid
        if ~validPSON(str)
            str = '';
            fprintf(['****************************\n']);
            fprintf(['Scanner didn''t get a full read.\n']);
            fprintf(['Please scan again.\n']);
            fprintf(['****************************\n']);
        end
    end
end