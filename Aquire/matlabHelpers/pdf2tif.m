function [newFileName,I] = pdf2tif(fileName)

    newFileName = strrep(fileName,'.pdf','.png');
    cmd = ['convert -density 150 -antialias "' fileName '" -append -resize 1024x -quality 100 "' newFileName '"'];
    system(cmd);
    if nargout == 2;I = imread(newFileName);end
    
end