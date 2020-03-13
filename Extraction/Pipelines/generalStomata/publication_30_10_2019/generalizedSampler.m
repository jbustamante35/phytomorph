function [subI] = generalizedSampler(I,T,D,SZ,disp,figH)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % I:= image to sample - one image
    % T:= transformation - one tranformation
    % D:= domains in cell array - multiple (cell array) domains
    % SZ:= domain size for reshape
    % disp:= for display mode
    % figH:=figure for display
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this will sample an image I with tranformation T over domain(s) D
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % assume froze tensor - this will break at some point and
    % it will make me unhappy - i should laugh - 31.10.2019
    if size(I,2) == 1
        I = thawTensor(I);
    end
    
    if ischar(I)
        I = generalizeLoader(I);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin < 5;disp = false;end
    if ((nargin < 6) && disp);figH = figure;end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % preallocate cell array
    subI = cell(numel(D)+1,1);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if T is displacement only
    % make it affine
    if numel(T) == 2
        dX = T(:)';
        T = eye(3);
        T(1:2,3) = dX;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make D affine if needed
    for d = 1:numel(D)
        if size(D{d},2) ~= size(T,2)
            D{d} = [D{d} ones(size(D{d},1),1)];
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % stack transformation/location for sampling
    subI{1} = T;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for each domains
    for d = 1:numel(D)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % transform D to tmpD
        tmpD = (T*D{d}')';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % clear tmpF
        tmpF = zeros([size(tmpD,1) size(T,3)]);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % for each channel
        for k = 1:size(I,3)
            % interpolation
            tmpF(:,k) = ba_interp2(I,tmpD(:,2),tmpD(:,1));
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if run with display on
        if disp
            figure(figH);
            imshow(I,[]);
            hold on
            plot(tmpD(:,2),tmpD(:,1),'.');
            hold off
            drawnow
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % reshape sample patches
        tmpF = reshape(tmpF,[SZ{d} size(tmpF,2)]);
        % store tensor for freeeze and transport
        subI{1+d} = tmpF;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % stack and freeze tensor
    subI = freezeTensor(subI);
end





