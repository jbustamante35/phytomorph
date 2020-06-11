function [TM,iTM] = ndmap(n,l)

    if nargin == 1;l=0;end
    
    % make hypersphere
    %%%%%%%%%%%%%%%%%%%%%
    [f{n-1},TM{n-1},x] = hSphere(n);
    x = sym(['x' num2str(n) '_'],[1 (n-1)]);
    for e = 1:numel(x);xn{e}=x(e);end
    
    % re-parameterize
    %%%%%%%%%%%%%%%%%%%%%
    TM{n-1} = symfun(TM{n-1}(1,xn{:}),x);
    
    
    % if there is a basement - go down the stairs
    %%%%%%%%%%%%%%%%%%%%%
    if n >= 3
        % call next level down
        %%%%%%%%%%%%%%%%%%%%%
        [tmp] = ndmap(n-1,l+1);
        for e = 1:numel(tmp)
            % get the map from below 
            TM{e} = tmp{e};
        end
    end
    
    
    % redo vars
    %%%%%%%%%%%%%%%%%%%%%
    X = [];
    for e = 1:numel(TM)
        X = [X,symvar(TM{e},'x')];
    end
    X = unique(X);
    
    % block to upscale dim
    %%%%%%%%%%%%%%%%%%%%%
    a = zeros(n,1);
    b = zeros(1,n+1);b(1) = 1;
    for e = 1:numel(TM)
        TM{e} = symfun([b;[a,TM{e}]],X);
    end
    
    if l == 0
        rot = TM{1};
        for e = 2:numel(TM)
            rot = TM{e}*rot;
        end
        
        % make variables
        %%%%%%%%%%%%%%%%%%%%%
        d = sym('d',[1 n]);
        s = sym('s',[1 n]);
        w = [X d s];
        % make T(w)
        %%%%%%%%%%%%%%%%%%%%%
        rot = symfun(rot,w);
        % make displacement
        %%%%%%%%%%%%%%%%%%%%%
        dx = [];
        for e = 1:n
            dx = [dx ; symfun(d(e),w)];
        end
        dx = [[1 zeros(1,n)];[dx symfun(eye(n),w)]];
        % make strain
        %%%%%%%%%%%%%%%%%%%%%
        z = zeros(1,n+1);
        z(1) = 1;
        strain = [z;[zeros(n,1),symfun(diag(s),w)]];
            
        
        total = dx*rot*strain;
        
        
        forwardMap = matlabFunction(total);
        backwardMap = matlabFunction(inv(total));
        TM = vectorizeFunc(forwardMap);
        iTM = vectorizeFunc(backwardMap);
    end
    
    %{
    % make displacement
    %%%%%%%%%%%%%%%%%%%%%
    dx = [];
    for e = 1:n
        dx = [dx ; symfun(d(e),[x s d])];
    end
    %}
    
    
    %{
    % make variables
    %%%%%%%%%%%%%%%%%%%%%
    d = sym('d',[1 n]);
    s = sym('s',[1 n]);
    % make displacement
    %%%%%%%%%%%%%%%%%%%%%
    dx = [];
    for e = 1:n
        dx = [dx ; symfun(d(e),[x s d])];
    end
    % make strain
    %%%%%%%%%%%%%%%%%%%%%
    z = zeros(1,n+1);
    z(end) = 1;
    S = [[symfun(diag(s),[x s d]),zeros(n,1)];z];
    TM = [[symfun(TM,[x s d]),dx];[zeros(1,nd) 1]];
    TM = TM*S;
    
    gTM = matlabFunction(inv(TM));
    gTM = vectorizeFunc(gTM);

    forwardMap = vectorizeFunc(matlabFunction(TM));
    backwardMap = vectorizeFunc(matlabFunction(inv(TM)));
    %}
end