function [dP] = extractLevelsFromContourMatrix(C)

        str = 1;
        cnt = 1;
        dP = {};
        while str < size(C,2)
            stp = str + C(2,str);
            dP{cnt} = C(:,(str+1):stp);
            cnt = cnt + 1;
            str = stp + 1;
        end
        
end