function [Y] = rpy(theta)
    CL = {'r','r','r'};
    Y = eye(3);
    X = eye(3);
    for k = 1:3
        quiver3(0,0,0,Y(1,k),Y(2,k),Y(3,k),CL{k});hold on
    end
   
    Y = rotW(X,Y,[1 2],theta(1));
    CL = {'g','g','g'};
    for k = 1:3
        quiver3(0,0,0,Y(1,k),Y(2,k),Y(3,k),CL{k});hold on
    end
    
    Y = rotW(X,Y,[2 3],theta(2));
    
    CL = {'b','b','b'};
    for k = 1:3
        quiver3(0,0,0,Y(1,k),Y(2,k),Y(3,k),CL{k});hold on
    end
    
    Y = rotW(X,Y,[1 2],theta(3));
    
    CL = {'k','k','k'};
    for k = 1:3
        quiver3(0,0,0,Y(1,k),Y(2,k),Y(3,k),CL{k});hold on
    end
end


function [Y] = rotW(X,Y,planeN,angle)

    rotM = eye(size(X));
    M = [[cos(angle),-sin(angle)];[sin(angle),cos(angle)]];
    rotM(planeN(1),planeN(1)) = M(1,1);
    rotM(planeN(1),planeN(2)) = M(1,2);
    rotM(planeN(2),planeN(1)) = M(2,1);
    rotM(planeN(2),planeN(2)) = M(2,2);
    Y = X*rotM*X'*Y;
end

%{
    close all
    X = eye(3);
    Y = rpy([pi/8 pi/16 pi/4]);
%}