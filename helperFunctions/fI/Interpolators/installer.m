function installer(installfolder)

    cd(installfolder);
    % make object file
    mex interpUtil.cpp -c
    % link
    mex ppmval.cpp interpUtil.obj 
    % compiling ppuval
    mex ppuval.cpp interpUtil.obj
    
end