function [Res] = getScreenResolution()
    %Sets the units of your root object (screen) to pixels
    set(0,'units','pixels') ;
    %Obtains this pixel information
    Pix_SS = get(0,'screensize');
    %Sets the units of your root object (screen) to inches
    set(0,'units','inches');
    %Obtains this inch information
    Inch_SS = get(0,'screensize');
    %Calculates the resolution (pixels per inch)
    Res = Pix_SS./Inch_SS;
    Res = Res(3:4);
end