function [FileList] = locateNumberedImageSets(ext)
    if nargin == 0;ext = {'TIF','JPG','JPEG','JPEG','TIF','BMP','PNG'};end
    % build extension string
    ext = cat(2,ext,lower(ext));
    EXT = '(';
    for e = 1:numel(ext)
        EXT = [EXT ext{e} '|'];
    end
    EXT(end) = ')';
    tic;
    CMD = ['locate -d /mnt/scratch4/spaldingdata.db -d /mnt/scratch4/tetra.db -d /mnt/scratch4/spaldingimages.db -d /mnt/scratch4/snapper.db --regex ''\/[0-9]*\.' EXT '$'''];
    [r,o] = system(CMD);
    FileList = splitlines(o);
    fprintf(['Located:' num2str(numel(FileList)) ' images in:' num2str(toc) '.\n']);
end