function [data_out] = generateAffineSequence(data_in,p,binaryDomain)
    % data_in(1,:)  := raw signal trace
    % data_in(2,:)  := normalized gradient_d1
    % data_in(2,:)  := normalized gradient_d2
    
    tidx = find(binaryDomain);
    for e = 1:numel(tidx)
        data_out(:,1,e) = [data_in(2:3,tidx(e));0;0];
        data_out(:,2,e) = [data_out(2,1,e);-data_out(1,1,e);0;0];
        data_out(:,3,e) = [0;0;0;0];
        data_out(:,4,e) = [p(2);p(1);tidx(e);1];
    end
end
    