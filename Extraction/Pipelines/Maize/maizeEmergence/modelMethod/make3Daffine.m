function [T] = make3Daffine(v)
    scale = v(1:3);
    rot = v(4:6);
    dx = v(7:9);
    
    S = diag([scale 1]);
    
    R1 = [1 0 0 0; 0 cos(rot(1)) -sin(rot(1)) 0; 0 sin(rot(1)) cos(rot(1)) 0; 0 0 0 1];
    R2 = [cos(rot(2)) 0 sin(rot(2)) 0; 0 1 0 0; -sin(rot(2)) 0 cos(rot(2)) 0;0 0 0 1];
    R3 = [cos(rot(3)) -sin(rot(3)) 0 0; sin(rot(3)) cos(rot(3)) 0 0; 0 0 1 0;0 0 0 1];
   
    D = eye(4);
    D(1:3,4) = dx(:);
    
    T = D*S*R1*R2*R3;
end

%{
    [n1,n2,n3] = ...
    ndgrid(linspace(-5,5,10),linspace(-5,5,10),linspace(-5,5,10));
    X = [n1(:) n2(:) n3(:) ones(size(n1(:)))];
    initv = [1 1 1 0 0 0 0 0 0 ];
  
    % scale test
    ds = linspace(1,5,5);
    close all
    for testEle = 1:3
        for e = 1:numel(ds)
            v = initv;
            v(testEle) = ds(e);
            a = make3Daffine(v);
            nX = (a*X')';
            plot3(nX(:,1),nX(:,2),nX(:,3),'.')
            axis([-10 10 -10 10 -10 10]);
            drawnow
            pause(.5);
        end
    end
    % rot test
    close all
    ds = linspace(-pi,pi,20);
    for testEle = 1:3
        for e = 1:numel(ds)
            v = initv;
            v(testEle+3) = ds(e);
            a = make3Daffine(v);
            nX = (a*X')';
            plot3(nX(:,1),nX(:,2),nX(:,3),'.')
            axis([-10 10 -10 10 -10 10]);
            view([0 90]);
            drawnow
            pause(.5);
        end
    end
    % displacement test
    close all
    ds = linspace(-2,2,10);
    for testEle = 1:3
        for e = 1:numel(ds)
            v = initv;
            v(testEle+6) = ds(e);
            a = make3Daffine(v);
            nX = (a*X')';
            plot3(nX(:,1),nX(:,2),nX(:,3),'.')
            axis([-10 10 -10 10 -10 10]);
            view([0 90]);
            drawnow
            pause(.5);
        end
    end

%}