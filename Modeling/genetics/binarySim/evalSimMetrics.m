function [d] = evalSimMetrics(programOutPut,toMatchY,type)
    switch type
        case 'hamming'
            d = bitwise_hamming(programOutPut,toMatchY);
           
        case 'lda'
            %UQ = unique(programOutPut);
            UQ = [0 1];
            for u = 1:numel(UQ)
                idx = find(programOutPut == UQ(u));
                data = toMatchY(idx);
                U(u) = mean(data);
                S(u) = std(data)^2;
            end
            d = -abs(diff(U))/mean(S);
            if isnan(d)
                d = 0;
            end
        case 'information'
             d = -mutualinfo(programOutPut,toMatchY);
        case 'matt'
             d = -matthews_correlation(toMatchY,programOutPut);
            
            
    end
end