function [out] = daqThread(channelToClient)
    
    backChannel = parallel.pool.PollableDataQueue;

    %{
    q = qConstant.Value;
    flag = true;
    while flag
        msg = poll(q,Inf);
        if isa(msg,'cell')
        if strcmp(msg,'Stop')
            flag = false;
        end
    end
    %}
    %out = labindex;
    %pause(30);
    %out = getCurrentWorker();
    %out = getLocalPart(qConstant);
    send(channelToClient,backChannel);
    out = 'hello';
    flag = true;
    while flag
        data = poll(backChannel,inf);
        if strcmp(data,'stop')
            flag = false;
        end
    end
end