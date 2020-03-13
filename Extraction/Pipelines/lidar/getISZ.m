function [iSZ,offset] = getISZ(H,W,TM,mag)


    rpH = round(mag*H);
    rpW = round(mag*W);
    rpT = round(mag*TM);
    
    offset = [min(rpH(:)) - 1 min(rpW(:)) - 1 min(rpT(:)) - 1];
    

    rpH = rpH - offset(1);
    rpW = rpW - offset(2);
    rpT = rpT - offset(3);
    
    iSZ = [max(rpH(:)),max(rpT(:))];
end