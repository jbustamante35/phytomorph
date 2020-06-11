% basis T 
% this is the basic block for building out tensors
%% load project
project.load('basisT');
%% save project
project.save();
%% mod-date: April 17, 2020
% thinking about 1) basis tensors as a function of a) vectors/parameters b)
% other tensors as an extension 2) mtimes as collapse of index notation - i
% will attempt to implement this as A*(B,2)
% Ended up useing [3,1]*[A,B];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test multiplication
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = basisT(rand(6,7,4));
B = basisT(rand(3,2,7));
C = [2;3]*[A,B];
%%
A = basisT(rand(6,7,4));
B = basisT(rand(4,2,7));
C = [[2;3],[3;1]]*[A,B];
%% generate a tensor of tensors
clear T
n = 20;
for i = 1:4
    for j = 1:5
        for k = 1:20
            d(i,j,k,:,:) = rand(n);
            T(i,j,k) = basisT(squeeze(d(i,j,k,:,:)));
        end
    end
end
sum(T,1)

