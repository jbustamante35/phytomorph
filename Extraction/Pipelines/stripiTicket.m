function [rawName,ticket] = stripiTicket(fileName,dispMode)
    if nargin == 1
        dispMode = true;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % strip ticket data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if dispMode;fprintf(['***********************************\n']);end
    if dispMode;fprintf(['Start ticket strip:' fileName '\n']);end
    if dispMode;fprintf(['***********************************\n']);end
    % if the code was passed a file name
    if ~isempty(fileName)
        % added to handle # in the file name
        [pth,nm,ext] = fileparts(fileName);
        % find the # symbol
        fidx = strfind(ext,'#');
        if ~isempty(fidx)
            rawName = [pth nm ext(1:(fidx(1)-1))];
            ticket = ext((fidx(1)+1):(fidx(2)-1));
            % remove and replaced with the above on 10.30.2019
            %rawName = fileName(1:(fidx(1)-1));
            %ticket = fileName((fidx(1)+1):(fidx(2)-1));
        else
            % if the extension does not have #
            rawName = fileName;
            ticket = '';
        end
    else
        rawName = fileName;
        ticket = '';
    end
    if dispMode;fprintf(['***********************************\n']);end
    if dispMode;fprintf(['End ticket strip:' fileName '\n']);end
    if dispMode;fprintf(['Raw Name-->' rawName '\n']);end
    if dispMode;fprintf(['Ticket Data-->' ticket '\n']);end
    if dispMode;fprintf(['***********************************\n']);end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % strip ticket data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end