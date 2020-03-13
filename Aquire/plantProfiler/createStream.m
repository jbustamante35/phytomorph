function [] = createStream(program,type,streamID,target,source)
    CMD = ['createStream ' program ' ' type ' ' num2str(streamID) ' ' num2str(1) ' ' target];
    if nargin == 5;CMD = [CMD ' ' source];end
    system(CMD,'-echo');
end