function [stateSpark] = sparkProcess(commandString,stateSpark,stripState)
    stateSpark = addToStruct(commandString,stateSpark);
end