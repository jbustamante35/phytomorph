function [newName] = recycleObject(object,arrayField)
    if nargin == 1;arrayField=[];end
    
    
    if ~isempty(arrayField)
        recycleObjectArray(object.(arrayField));
    end
    
    objectHome = getenv('PHYTO_OBJECTS');
    recyclePath = [objectHome '.recycle' filesep];
    mkdir(recyclePath);
    
    % get the file if it is not already
    if ~isa(object,'file')
        % if pointer then we want to remove the object it points to
        % not the pointer
        if isPtr(object)
            object = object.refs;
        end
        % get file name
        [fileName] = getObjectFile(object);
    end
    
    
    [~,nm] = fileparts(fileName);
    newName = [recyclePath nm];
    %movefile(fileName,newName);
end