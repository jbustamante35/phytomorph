% make e - element is size [3 3];
% array of size [1 1 4 9];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear e
for i = 1:4
    for j = 1:9
        P = rand(1,3);
        E = rand(3,3);
        e(1,1,i,j) = basisT(E,P,[]);
    end
end
size(e)
%% make f - element is size [1 3 2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear f
for i = 1:3
    for j = 1:2
        P = rand(1,3);
        E = rand(1,3,2);
        f(i,j) = basisT(E,P,[]);
    end
end
size(f);
%% add a = e + f
% at the scaler level test this first
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e_rand = rand(size(e));
f_rand = rand(size(f));
a_rand = e_rand + f_rand;
sza_rand = size(a_rand);
e_vec = e_rand(:);
f_vec = f_rand(:);
a_vec = e_vec' + f_vec;
for i = 1:numel(a_vec)
    per(i) = find(a_vec(i) == a_rand(:));
end
%% this takes some time to run - let it be
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e_rand = e;
f_rand = f;
a_rand = e_rand + f_rand;
sza_rand = size(a_rand);
e_vec = e_rand(:);
f_vec = f_rand(:);
a_vec = e_vec' + f_vec;
for i = 1:numel(a_vec)
    per(i) = find(a_vec(i) == a_rand(:));
end
a_rand = e_rand + f_rand;
sza_rand = size(a_rand);
e_vec = e_rand(:);
f_vec = f_rand(:);
a_vec = e_vec' + f_vec;
for i = 1:numel(a_vec)
    perT(i) = find(a_vec(i) == a_rand(:));
end
%% try recursion - this might reveal - dim,rank,level
% thnk about 2D array, matrix, 2nd order/way tensor
% the elements are 1st order/way/rank, vectors of dim 20
% then A(i,j).E(k) could be A(i,j,k) if A(i,j) is homogenous

for l = 1:5
    for k = 1:9
        clear f
        for i = 1:3
            for j = 1:2
                P = rand(1,3);
                E = rand(1,3,2);
                f(i,j) = basisT(P,[],E);
            end
        end
        P = [];
        OUT(l,k) = basisT(P,[],f);
    end
end


%%
clear g
for i = 1:3
    for j = 1:2
        P = rand(1,3);
        E = rand(3,1);
        g(i,j) = basisT(P,[],E);
    end
end
clear a b x
a = e.applyOperation(@(x,y)plus(x,y),f);
b = e.applyOperation(@(x,y)mtimes(x,y),g);
e.applyOperation(@(x)x(:));
%%

