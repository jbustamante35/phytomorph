function [fld] = extractFieldFromPath(fileName,N)
    [p,n,ext] = fileparts(fileName);
    p = [p filesep];
    pidx = strfind(p,filesep);
    str = pidx(end-N)+1;
    stp = pidx(end-(N-1))-1;
    fld = p(str:stp);
end