function msg = hdferrmsg(status,interface,function_name)
%Use the last error generated by the HDF library to generate an error message

%   Copyright 2007-2011 The MathWorks, Inc.
%   $Revision: 1.1.8.4 $  $Date: 2011/05/17 02:27:22 $

if status==-1
    error_code = hdfhe('value',1);
    if strncmpi(hdfhe('string',error_code),'no error',8)
        error(message('MATLAB:imagesci:hdf:unknownLibraryError', function_name, interface));
    else
        error(message('MATLAB:imagesci:hdf:libraryError', hdfhe( 'string', error_code )));
    end
end
