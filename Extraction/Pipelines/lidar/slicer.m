function [I,offset] = slicer(H,W,TM,WTH,SZ,offset,mag)

    pidx = W > WTH(1) & W < WTH(2);
    pH = H(pidx);
    pW = W(pidx);
    pTM = TM(pidx);

    %%
    rpH = round(mag*pH);
    rpW = round(mag*pW);
    rpT = round(mag*pTM);
   
    
    rpH = rpH - offset(1);
    rpW = rpW - offset(2);
    rpT = rpT - offset(3);

    I = zeros(SZ);
    IDX = sub2ind(size(I),rpH,rpT);
    numel(unique(IDX))
    numel(IDX)
    
    I(IDX) = rpW;

end