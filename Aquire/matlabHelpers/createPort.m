function [] = createPort(portID)
     CMD = ['createPort ' portID];
     system(CMD);
end