function tf = isras(filename)
%ISRAS Returns true for a RAS file.
%   TF = ISRAS(FILENAME)

%   Copyright 1984-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2005/11/15 01:12:41 $

fid = fopen(filename, 'r', 'ieee-be' );
if (fid < 0)
    tf = false;
else
    sig = fread(fid, 1, 'uint32');
    fclose(fid);
    tf = isequal(sig, 1504078485);      % 0x59a66a95
end
