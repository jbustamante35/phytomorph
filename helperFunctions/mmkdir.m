function [] = mmkdir(p)
    CMD = ['mkdir -p ' p];
    system(CMD);
end