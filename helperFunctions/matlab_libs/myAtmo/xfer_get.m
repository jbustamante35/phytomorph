function [localList] = xfer_get(remoteList,localPath,disjoint_union,disp,ticketList)
    try
        %%%%%%%%%%%%%%%%%%
        % remote list: - cell array of strings for files
        % local path: - string pointing to local location to put files
        % disjoint_union: - if true-place prefix of number in local copy
        % disp: if true - then display
        %%%%%%%%%%%%%%%%%%
        % set is cel for return
        isCell = isa(remoteList,'cell');
        % assume tmp local path if not provided
        if nargin == 1;localPath = makeTempLocation();end
        % assume tmp local path if localPath is empty
        if isempty(localPath);localPath = makeTempLocation();end
        % set disp
        if nargin < 4;disp = false;end
        % demand cell input
        if ~iscell(remoteList);remoteList = {remoteList};end
        
        
        
        %%%%%%%%%%%%%%%%%%
        % xfer from irods OR idfunction for local
        % uses icommands to tranfer to local file system
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if disjoint union is null
        % default to false
        if nargin < 3;disjoint_union = false;end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if ticket is given - give location for ticket
        if nargin <= 4
            ticketCMD = '';
        else
            ticketCMD = ' -t "<T>" ';
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % create template commands
        % hardcoded force overwrite
        cmd_template = ['iget -f -N 2 -V' ticketCMD ' "<S>" "<D>"'];
        %cmd_template = ['iget -f -V' ticketCMD ' "<S>" "<D>"'];
        local_template = [localPath '<D>'];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make the target directory
        mkdir(['.' localPath]);
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i = 1:numel(remoteList)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if it is a iRODS file then transfer to localPath
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isIRODS(remoteList{i})
                
                
                [fileName,ticket] = stripiTicket(remoteList{i});
        
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % report with wait bar
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if disp
                    if i == 1
                        h = waitbar(i/numel(remoteList),'Transferring data from iRODS to local machine');
                    else
                        waitbar(i/numel(remoteList),h,'Transferring data from iRODS to local machine');
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % xfer file
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [pth,nm,ext] = fileparts(fileName);
                cmd = strrep(cmd_template,'<S>',fileName);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % if disjoin union flag - then tag with i-number
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if disjoint_union
                    cmd = strrep(cmd,'<D>',[localPath num2str(i) '-' nm ext]);    
                else
                    cmd = strrep(cmd,'<D>',[localPath nm ext]);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % if ticketCMD isempty
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if ~isempty(ticket)
                    cmd = strrep(cmd,'iget ',['iget -t ' ticket ' ']);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % perform transfer
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if ~isdeployed
                    [status, result] = system(cmd);
                end
                
                if isdeployed
                    fprintf(['System call:' cmd '\n']);
                    [status, result] = system(cmd,'-echo');
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % report status of operation
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if status == 0
                    fprintf(['xFer:success: ' [nm ext] '\n']);
                else
                    fprintf(['xFer:fail: ' [nm ext] '\n']);
                    fprintf(['command@' cmd '\n']);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % get local file name for return
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if disjoint_union
                    localList{i} = strrep(local_template,'<D>',[num2str(i) '-' nm ext]);
                else
                    localList{i} = strrep(local_template,'<D>',[nm ext]);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if it is a localfile then do nothing
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif isWID(remoteList{i})
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % strip down file
                [pth,nm,ext] = fileparts(remoteList{i});
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % report with wait bar
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if disp
                    if i == 1
                        h = waitbar(i/numel(remoteList),'Transferring data from CHTC to local machine');
                    else
                        waitbar(i/numel(remoteList),h,'Transferring data from CHTC to local machine');
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % if disjoin union flag - then tag with i-number
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if disjoint_union
                    localName = [localPath num2str(i) '-' nm ext];  
                else
                    localName = [localPath nm ext];
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                localList{i} = localName;
                % create copy command
                cmd = ['mc cp ' remoteList{i}(2:end) ' ' localName]; 
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % perform transfer
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if ~isdeployed
                    [status, result] = system(cmd);
                end
                
                if isdeployed
                    fprintf(['System call:' cmd '\n']);
                    [status, result] = system(cmd,'-echo');
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % report status of operation
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if status == 0
                    fprintf(['xFer:success: ' [nm ext] '\n']);
                else
                    fprintf(['xFer:fail: ' [nm ext] '\n']);
                    fprintf(['command@' cmd '\n']);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % is local file - pass through
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                localList{i} = remoteList{i};
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % unwrap - maybe
        if (~(isCell) && (numel(localList) == 1))
            localList = localList{1};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % close wait bar
        if disp;close(h);end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % return the local list
        if iscell(localList);locallist = localList';end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    catch ME
        errorReport = getReport(ME);
        errorReport
        fprintf(['error@xfer_get\n']);
    end
end