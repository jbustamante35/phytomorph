function [F,K] = smoothTasselContour(J,para)
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    About:      
                cwtK_closed_imfilter.m  (Inputs are relative to 1200dpi)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dependency: 
                fspecial.m, gradient.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Variable Definition:
                J:      The information is needed. 
                para:      The information is needed.
                
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    % MUST REMOVE THE BASE
    % FIND AND SUPPRES    
    try
        windowMAG= 5;
        h0 = fspecial('gaussian',[windowMAG*para 1],para);
        h0 = h0 / sum(h0);
        % calculate curvature            
        tmp = imfilter(J,h0,'circular');
        d1X = gradient(tmp')';
        d2X = gradient(gradient(tmp'))';            
        K = (d1X(:,1).*d2X(:,2) - d1X(:,2).*d2X(:,1)).*(d1X(:,1).^2 + d1X(:,2).^2).^-3/2;    
       
        for e = 1:size(d1X,1)
            T = d1X(e,:);
            T = T / norm(T);
            N = d2X(e,:);
            N = N / norm(N);
            %quiver(J(e,1),J(e,2),T(1),T(2),'g')
            F(:,:,e) = [T;N];
        end
        F = permute(F,[3 2 1]);
        
    catch ME
        close all;
        getReport(ME);
        fprintf(['******error in:cwtK_closed_imfilter.m******\n']);
    end

end