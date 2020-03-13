function [stop] = outputPLT(func,A,B,C)
    if isstruct(A)
        pt = A.bestx;
    else
        pt = A;
    end
    
    if isstruct(A)
        iter = A.iteration;
    else
        iter = B.iteration;
    end
    [~,bd] = func(pt);
    global BD
    BD = [BD;bd];
    
    if iter > 10
        CL = {'r' 'g' 'b'};
        for e = 1:size(BD,2)
            plot(BD(10:end,e),CL{e})
        end
    end
    hold on
    stop = false;
end