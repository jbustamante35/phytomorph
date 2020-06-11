function [hs,nTM,x] = hSphere(n,hs,l)
    %%%%%%%%%%%%%%%%%%%%%%%%
    % return:
    % 1) function for position on hypersphere
    % 2) function for normal space + tangent space
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % recursive build up for hypersphere
    %%%%%%%%%%%%%%%%%%%%%%%%
    if nargin == 1
        l = 1;
        x = sym('x',[1 (n)]);
        hs = symfun(x(1),x);
        hs = hSphere(n,hs,l);
    else
        if l ~= n
            x = symvar(hs,'x');
            z = symfun([zeros(1,l-1),sin(x(l+1))],x);
            e = eye(l);
            M = symfun(reshape([e(1:(end-1))';cos(x(l+1))],[l,l]),x);
            M = [M;z];
            hs = M*hs;
            l = l + 1;
            hs = hSphere(n,hs,l);
        end
    end
    
    % build up for normal and tangent space
    %%%%%%%%%%%%%%%%%%%%%%%%
    % get input values
    x = symvar(hs,'x');
    % assume all are readl
    assume(x,'real');
    % assume radius is positive
    assume(x(1),'positive');
    % assume 2 -> (n-1) is [0,pi]
    for e = 2:(numel(x)-1)
        assume(x(e) >= 0 & x(e) <= pi);
    end
    assume(x(end) >= -pi & x(end) <= pi);

    % build out the basis vectors as column vectors
    TM = [];
    for e = 1:numel(x)
        TM = [TM diff(hs,x(e))];
    end
    TM = symfun(TM,x(2:(end)));
    nTM = [];
    for e = 1:numel(x)
        tmp = diff(hs,x(e));
        ntmp = sum((tmp.*tmp)).^.5;
        ntmp = simplify(ntmp,'Steps',3);
        tmp = tmp*ntmp^-1;
        tmp = simplify(tmp,'Steps',3);
        nTM = [nTM tmp];
    end
    nTM = symfun(nTM,x(2:(end)));
    
end
%{
    %%%%%%%%%%%%%%%%%%%%%
    nd = 3;
    [f,TM] = hSphere(nd);
    g = matlabFunction(f);
    gTM = matlabFunction(TM);





    close all
    r = 1;t1 = linspace(0,pi,30);t2 = linspace(-pi,pi,20);
    for i = 1:numel(t1)
        p = [];
        for e = 1:numel(t2)
            p(:,e) = g(r,t1(i),t2(e));
        end
        plot3(p(1,:),p(2,:),p(3,:),'k');hold on
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CL = {'r','g','b'};
    r = 1;t1 = 0;t2 = linspace(-pi,pi,20);
    for e = 1:numel(t2)
        p = g(r,t1,t2(e));
        T = gTM(r,t1,t2(e));
        for k = 1:size(T,2)
            quiver3(p(1),p(2),p(3),T(1,k),T(2,k),T(3,k),CL{k});
            drawnow
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CL = {'r','g','b'};
    r = 1;t1 = pi + pi/8;t2 = linspace(-pi,pi,20);
    for e = 1:numel(t2)
        p = g(r,t1,t2(e));
        T = gTM(r,t1,t2(e));
        for k = 1:size(T,2)
            quiver3(p(1),p(2),p(3),T(1,k),T(2,k),T(3,k),CL{k});
            drawnow
        end
    end


    [P] = randomPDF(3);

    
    

%}