function [fileName,iticket] = issueTicket(fileName,uses,type)
    % default to 1-use read tickets 
    if nargin == 1;uses = 1;end
    if nargin <= 2;type = 'read';end
    iticket = '';
    if isIRODS(fileName)
        % create the ticket
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cmd = ['iticket create ' type ' ' '"' fileName '"'];
        [o,r] = system(cmd);
        fidx = strfind(r,':');
        iticket = r(fidx+1:end-1);
        
        % modify uses
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if uses ~= Inf
            cmd = ['iticket mod ' iticket ' uses ' num2str(uses)];
            %[o,r] = system(cmd);
        end


        % make write ticket 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(type,'write')
            cmd = ['iticket mod ' iticket ' write-file ' num2str(uses)];
            cmd = ['iticket mod ' iticket ' write-file ' num2str(0)];
            [o,r] = system(cmd);
        end
        
        % add expire date
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        expireDate = datestr(addtodate(now,11,'day'),'YYYY-mm-dd.hh:MM:ss');
        cmd = ['iticket mod ' iticket ' expire ' expireDate];
        [o,r] = system(cmd);
        
        % return file name
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fileName = [fileName '#' iticket '#'];
    elseif isWID(fileName)
        fileName(1) = [];
        cmd = ['mc share download ' fileName];
        [o,result] = system(cmd);
        result = strtrim(result);
        key = 'Share: ';
        fidx = strfind(result,key);
        fileName = result((fidx(1)+numel(key)):end);
        
    end
end