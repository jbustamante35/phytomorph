function [] = mountCyVerse(source,target)
    uname = input('CyVerse User Name:','s');
    pw = input('CyVerse Password:','s');
    initStr = [uname '\n' pw '\n'];
    CMD = ['echo ' initStr ' | sudo ~/mountCyverse ' source ' ' target];
    [o,r] = system(CMD);
end