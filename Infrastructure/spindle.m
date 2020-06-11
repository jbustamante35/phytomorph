function [] = spindle(data,file)
    spoolNumeric(data,file)
end


function [] = spoolNumeric(data,file)

        

    type = class(data);
    
    
    tid = fopen(file,'w');
    
    
    % write the escape char
    escapeMax = 2^8;
    escapeType = 'uint8';
    escape = randi(escapeMax,escapeType);
    %escape = typecast(escape,escapeType);
    fwrite(tid,escape,type,'n');
    %
    stopMax = 2^8;
    stopType = 'uint8';
    stop = randi(stopMax,stopType);
    %stop = typecast(stop,stopType);
    fwrite(tid,stop,type,'n');
    
    
    data = typecast(data,'uint64');
    % type cast to search for stop codes
    % and apply escape
    
    
    
    switch type
        case 'double'
            data = typecast(data(:),'uint8');
    end
    
    eidx = find(data == escape);
    sidx = find(data == stop);
    tidx = [eidx;sidx];
    
    newData = zeros([1,numel(data)+numel(eidx)+numel(sidx)],class(data));
    strX = 1;
    strY = 1;
    dtidx = diff(tidx);
    for e = 1:numel(tidx)
        
    end
    
    
    fwrite(tid,data,type,'n');
    fclose(tid);
end

function [] = spoolChar()
end

function [] = spoolCell()
end

%{
    N = randi(2000,1);
    data = rand(N,1);
    file = '~/data';
    spindle(data,file);
%}