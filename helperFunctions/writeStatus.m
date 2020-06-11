function [] = writeStatus(oPath,nm,status)
    if status == 1
        sFile = [oPath '{fileName_' nm '}{status_success}.txt'];
        status = num2str(status);
    else
        sFile = [oPath '{fileName_' nm '}{status_fail}.txt'];
        status = num2str(status);
    end
    fileID = fopen(sFile,'w');
    fprintf(fileID,'%s',status);
    fclose(fileID);
end