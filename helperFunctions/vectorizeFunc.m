function [func] = vectorizeFunc(func)
    sfunc = func2str(func);
    fidx1 = strfind(sfunc,'(');
    fidx2 = strfind(sfunc,')');
    args = sfunc((fidx1(1)+1):(fidx2(1)-1));
    body = sfunc((fidx2(1)+1):end);
    args = [',' args ','];
    fidx = strfind(args,',');
    for e = 1:(numel(fidx)-1)
        source{e} = args((fidx(e)+1):(fidx(e+1)-1));
        target{e} = ['x(' num2str(e) ')'];
    end
    for e = 1:numel(target)
        body = strrep(body,source{e},target{e});
    end
    sfunc = ['@(x)' body];
    func = str2func(sfunc);
end