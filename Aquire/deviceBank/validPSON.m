function [v] = validPSON(inString)
    if ~isempty(inString)
        for e = 1:numel(inString)
            str(e) = -strcmp(inString(e),'{');
            stp(e) = strcmp(inString(e),'}');
        end
        v = sum(str + stp) == 0;
    else
        v = false;
    end
end