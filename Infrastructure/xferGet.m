function [r] = xferGet(source,target,varargin)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the argumet list
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kvp = pargin(varargin{:});
    
    % if clip:source is given
    if any(contains({kvp.key},'clip:source'))
        
        % get the valueof clip:source
        kidx = find(strcmp({kvp.key},'clip:source'));
        % get the clip value
        clipValue = kvp(kidx).value;
        
        % if char
        if ischar(clipValue)
            % if *target then use the target to clip
            if strcmp(clipValue,'*target')
                clipValue = target;
            end
            % if '' then use '/' to make no clip
            if strcmp(clipValue,'')
                clipValue = filesep;
            end
            clipValue = numel(strfind(clipValue)) - 1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % build out the target values
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for e = 1:numel(source)
            [pth,~,~] = fileparts(source{e});
            sidx = strfind(pth,filesep);
            destinationPath = pth((sidx(clipValue+1)+1):end);
            targetList{e} = [target destinationPath filesep];
        end
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % replicate the target
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make target
        mmkdir(target);
        % replicate the target
        if ~iscell(target)
            for e = 1:numel(source)
                targetList{e} = target;
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find the unique targets
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    UQ = unique(targetList);
    for u = 1:numel(UQ)
        % make each target
        mmkdir(UQ{u});
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % execute commands
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for e = 1:numel(source)
        r(e) = subfXfer(source{e},targetList{e},varargin);
    end

end
