function [T,f] = findKernelPeriod(sig,N,MAXT,MXsupport)
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    About:      
                findT.m is main function to handle ear analysis. It takes all input variables 
                for its dependent functions. This function returns final result including image with 
                bounding box. (Inputs are relative to 1200dpi)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dependency: 
                bindVec.m, interp1.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Variable Definition:
                sig:       	
                N:           
                MAXT:          Input arguments passed.
                
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    T = NaN;
    f = NaN;
    try
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % high frequency cut off specifed as kernel period
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Tline = N./((1:numel(sig))-1)';
        T_thresh = Tline < MAXT;
        
        
        % changed on aug 11 2016 from 8 to 10
        [localMAX] = nonmaxsuppts(sig,MXsupport);

        % added for simulation
        %sig = sig.*T_thresh;

        %{
        nsig = bindVec(sig);
        thresh = graythresh(nsig);
        bidx = (nsig > thresh);
        %}
        %fidx = find(localMAX.*bidx.*T_thresh);
        
        fidx = find(localMAX.*T_thresh);
        
        fpeak = sig(fidx);
        [fpeak,sidx] = sort(fpeak,'descend');

        f = mean(fidx(sidx(1)));
        % first peak
        f = fidx(1);
        %%
        T = N/(f-1);
        f = T^-1;
    catch ME
        close all;
        getReport(ME);
        fprintf(['******error in:findT.m******\n']);
    end
end