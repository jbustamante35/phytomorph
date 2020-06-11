fileName = '/iplant/home/hirsc213/maizeData/seedlingData/10-May-2017/{Plot_894}{Experiment_44}{Planted_4-24-2017}{SeedSource_DI 2114-3}{SeedYear_2016}{Genotype_Mo17}{Treatment_Control}{PictureDay_16}.nef';
I = imread(fileName);
%%
G = rgb2gray(I);
%%
close all
b = bitget(G,8);
imshow(b,[]);
%%
sub = double(b(1:1000,1:1000));
c = im2colF(sub,[51 51],[1 1]);
%%
c = logical(bitget(c,1));
%% 
v = rand(size(c,1),10);
v = logical(round(v));
tic
d = pdist2(ct,v','hamming');
toc
%%
ct = c';
%%
options = optimoptions('particleswarm','SwarmSize',100,'UseVectorized',true,'Display','iter');
func = @(x)mean(pdist2(bitget(round(x),1),ct(1:1000,:),'hamming'),2);
d = func(v');
x0 = particleswarm(func,size(ct,2),zeros(size(ct,2),1),ones(size(ct,2),1),options);
%%
it = 100;
sub = ct(1:1000,:);
% get mean - first guess
u = logical(round(mean(sub,1)));
% make z vector
z = rand(it,size(ct,2)*2);
% bit get on z
z = bitget(round(z),1);

d = dmm(z,sub,u);

options = optimoptions('particleswarm','SwarmSize',100,'UseVectorized',true,'Display','iter');
func = @(z)dmm(z,sub,u);
x0 = particleswarm(func,size(ct,2)*2,zeros(size(ct,2)*2,1),ones(size(ct,2)*2,1),options);
%% how related are the measurements
sub = c(:,1:1000);
u = logical(round(mean(sub,2)));
sub = carveOP(sub,u);
dc = pdist2(sub,sub,'hamming');
[v,d] = eigs(dc,2);




