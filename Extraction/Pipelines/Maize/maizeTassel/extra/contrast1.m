function [d] = contrast1(dT,generatorT,sourceP,sourceGlob,sourceMask,targetF,targetM,x,szX)
    
    dT = generatorT(dT);
    
    sourceP = sourceP.push(dT);
    
    
    sourceF = sourceP.getF(sourceGlob,x,szX);
    %targetF = targetP.getF(targetGlob,x,szX);
    
    sourceM = sourceP.getF(sourceMask,x,szX);
    %targetM = targetP.getF(targetMask,x,szX);
    %sourceM = sourceM > .5;
    
    
    T = sourceF(:) - targetF(:);
    T = mtimesx(T,'T',T);
    B = mtimesx(targetM(:),'T',sourceM(:));
    d = T*B^-1;
    %d = T.^.5/B;
    %d = d - sum(B);
end