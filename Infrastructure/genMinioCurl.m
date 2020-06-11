function [curlCmd] = genMinioCurl(type,bucket,location,expire)
    cmd = ['mc --json share  ' type ' --recursive -E ' expire ' chtc/' bucket location];
    [o,result] = system(cmd);
    result = jsondecode(result);
    curlCmd = result.share;
end

%{
    curlCmd = genMinioCurl('upload','abc','/next/level/','10m');
%}