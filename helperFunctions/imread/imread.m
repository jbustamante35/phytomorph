function [X, map, alpha] = imread(varargin)
    % add class type file
    % add ticket get inside xfer_get
    % add ext=nms
    % add ext=nef
    % add JSON doid in clumps
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ischar(varargin{1})
        fprintf(['start:image read@' varargin{1} '\n']);tic
    elseif isa(varargin{1},'file')
        varargin{1} = [varargin{1}.filePath filesep varargin{1}.fileName];
        fprintf(['start:image read@ pass-through.\n']);tic
    end
    if ~isnumeric(varargin{1})

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % xfer_get the file if needed
        varargin{1} = xfer_get(varargin{1});
        % get file parts
        [pth,nm,ext] = fileparts(varargin{1});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % this is if the image file is a NEF
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmpi(ext,'.nef')
            % need to convert to TIF
            fprintf(['Converting NEF -> TIF \n']);
            % call dcraw
            CMD = ['dcraw -a -T "' varargin{1} '"'];
            %CMD = ['dcraw -T "' filename '"'];
            %CMD = ['dcraw -w -T "' filename '"'];
            system(CMD);
            %{
            if toDelete
                fprintf(['Deleting NEF\n']);
                delete(filename);
            end
            %}

            % pass on the file name
            varargin{1} = strrep(varargin{1},ext,'.tiff');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if nms file type
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif strcmp(ext,'.nms')
            % use the NMS converter
            s = openNMS(filename);

            %{
            a1 = double(fillmissing(s.z,'linear',1));
            a2 = double(fillmissing(s.z,'linear',2));
            X = .5*(a1+a2)/255;
            %}

            % convert to double and divide
            X = double(s.r)/255;
            % auto run adaptive histogram equal
            X = adapthisteq(X,'NumTiles',[15 15]);
        else
            % 
            [X, map, alpha] = imread_org(varargin{:});
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        X = varargin{1};
    end
    
    if ischar(varargin{1})
        fprintf(['end:image read@' num2str(toc) '\n']);
    end
end
