function [] = genUpload(remotePoint,exp)
    if nargin == 1;exp = '1d';end
    cmd = ['mc --json share upload --recursive -E ' exp ' chtc/' lower(remotePoint)];
    [o,result] = system(cmd);
    result = jsondecode(result);
    curlCMD = strrep(result.share,'<NAME>',"$remoteFile");
    curlCMD = strrep(curlCMD,'<FILE>',"$localFile");
    cpDataFileName = '/mnt/scratch1/phytomorph_dev/Infrastructure/cpData_template';
    cpFileData = fileread(cpDataFileName);
    cpFileData = strrep(cpFileData,'<curlCMD>',curlCMD);
    
    fileID = fopen('exp.txt','w');
fprintf(fileID,'%6s %12s\n','x','exp(x)');
fprintf(fileID,'%6.2f %12.8f\n',A);
fclose(fileID);
end