function [] = writeObject(object)
    [fileName] = getObjectFile(object);
    object = jsonencode(object);
    fileID = fopen(fileName,'w');
    fprintf(fileID,'%s',object);
end