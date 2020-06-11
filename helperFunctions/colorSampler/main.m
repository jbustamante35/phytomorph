%% test #0
% test the hyperSphere
[HS2,TM2] = hSphere(2);
[HS3,TM3] = hSphere(3);
[HS4,TM4] = hSphere(4);
%%
% for 2D
close all
mFunc = matlabFunction(HS2);
hs2 =  vectorizeFunc(mFunc);
mFunc = matlabFunction(TM2);
tm2 =  vectorizeFunc(mFunc);

[RAD,TH] = ndgrid(linspace(1,1,1), linspace(-pi,pi,100));
TH = [RAD(:) TH(:)];
for e = 1:size(TH,1)
    P = hs2(TH(e,:));
    plot(P(1),P(2),'k.');
    hold on
    
    W = tm2(TH(e,2));
    quiver(P(1),P(2),W(1,1),W(2,1),'r');
    quiver(P(1),P(2),W(1,2),W(2,2),'g');
    
    drawnow
    pause(.1)
end

%% for 3D
close all
mFunc = matlabFunction(HS3);
hs3 =  vectorizeFunc(mFunc);
mFunc = matlabFunction(TM3);
tm3 =  vectorizeFunc(mFunc);

[RAD,TH1,TH2] = ndgrid(linspace(1,1,1), linspace(pi/8,pi/8,1),linspace(-pi,pi,100));

TH = [RAD(:) TH1(:) TH2(:)];
for e = 1:size(TH,1)
    
    P = hs3(TH(e,:));
    plot3(P(1),P(2),P(3),'k.');
    hold on
    
    W = tm3(TH(e,2:3));
    quiver3(P(1),P(2),P(3),W(1,1),W(2,1),W(3,1),'r');
    quiver3(P(1),P(2),P(3),W(1,2),W(2,2),W(3,2),'g');
    quiver3(P(1),P(2),P(3),W(1,3),W(2,3),W(3,3),'b');
    drawnow
    pause(.1)
end



[RAD,TH1,TH2] = ndgrid(linspace(1,1,1), linspace(0,pi,50),linspace(-pi,pi,50));
TH = [RAD(:) TH1(:) TH2(:)];
for e = 1:size(TH,1)
    
    P = hs3(TH(e,:));
    plot3(P(1),P(2),P(3),'k.');
    hold on
    
    W = tm3(TH(e,2:3));
    drawnow
    
    W'*W
end




%% test #1
% Test with Julian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
viewDims = 2;
FileList= {};
FileList{1} = '/mnt/spaldingdata/nate/octerineDataStore/testData/channelSampler/200416_SegmentScores_3879Patches.mat';
dataLoader{1} = @(x,ni)matLoader(x,'SCRS',1:viewDims);
dataLoader{2} = @(x,ni,data)(1:size(data,1));
d = dataLoader{1}(FileList{1},1);
N = 30;
Nfiles = 1;
Npoints = 10*size(d,1);
range = [];N_spaceVecs = {};
for e = 1:size(d,2)
    range(e,:) = [min(d(:,e)) max(d(:,e))];
    N_spaceVecs{e} = linspace(range(e,1),range(e,2),N);
end
density = channelSampler(FileList,Nfiles,Npoints,N_spaceVecs,dataLoader);
%% show the test #1 results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
plot(d(:,1),d(:,2),'.')
waitforbuttonpress
close all
image = density.density;
image = interp2(image,2);
imshow(image,[]);
waitforbuttonpress
%% if 3D
%{
close all
plot3(d(:,1),d(:,2),d(:,3),'.')
waitforbuttonpress
close all
for e = 1:size(density.density,3)
    image = density.density(:,:,e);
    image = interp2(image,2);
    imshow(image,[]);
    waitforbuttonpress
end
%}
%% discusssion about next steps
% point #1 - fit univariate distribution with gaussian
close all
pd = fitdist(d(:,1),'norm');
ksdensity(d(:,1));hold on;
xi = linspace(min(d(:,1)),max(d(:,1)),1000);
p = pd.pdf(xi);
plot(xi,p,'g');
% point #2 - fit multivariate distribution with gaussian
% sample statistics are a good fit
waitforbuttonpress;close all
u = mean(d(:,1:2));
sig = cov(d(:,1:2));
pdfF = @(x)mvnpdf(x,u,sig);
NP = 1000;
xi = linspace(min(d(:,1)),max(d(:,1)),NP);
yi = linspace(min(d(:,2)),max(d(:,2)),NP);
[x,y] = ndgrid(xi,yi);
f = pdfF([x(:),y(:)]);
f = reshape(f,size(x));
imshow(f,[]);hold on;
nd = [];
nd(:,1) = NP*bindVec(d(:,1));
nd(:,2) = NP*bindVec(d(:,2));
xd = linspace(0,NP,NP);
yd = linspace(0,NP,NP);
[xx,yy] = ndgrid(xd,yd);
plot(nd(:,2),nd(:,1),'r.');
% point #3 - use Gaussian Mixture Model with mixture = 1
ops = statset('Display','iter');
GMModel = fitgmdist(d(:,1:2),1,'Options',ops);
f = GMModel.pdf([x(:),y(:)]);
f = reshape(f,size(x));
contour(xd,yd,f,5,'r');
% point #3 - use Gaussian Mixture Model with mixture = 2
ops = statset('Display','iter');
GMModel = fitgmdist(d(:,1:2),2,'Options',ops);
f = GMModel.pdf([x(:),y(:)]);
f = reshape(f,size(x));
contour(xd,yd,f,5,'b');
for e = 1:2
    u = GMModel.mu(e,:);
    sig = GMModel.Sigma(:,:,e);
    f = mvnpdf([x(:),y(:)],u,sig);
    f = reshape(f,size(x));
    contour(xd,yd,f,5,'c');
end
% point #4 - Univariate distributions can fit different families
waitforbuttonpress
close all
toFit = bindVec(d(:,1))+eps;
pd = fitdist(toFit,'Weibull');
ksdensity(toFit);hold on;
xi = linspace(0,1,1000);
p = pd.pdf(xi);
plot(xi,p,'g');
% point #5 - what do we do to fit multivariate mixtures of different types
% there is no tool to perform fit of weibull-gaussian
% point #6 - return to the question - how does a distribution get fit?
phat = mle(toFit,'distribution','Weibull');
clear data
contrastFunc = @(data,p,type)-sum(log(pdf(type,data,p(1),p(2))));
contrastFunc(toFit,phat,'Weibull')
func = @(para)contrastFunc(toFit,para,'Weibull');
xSol = fmincon(func,[1 1],[],[],[],[],zeros(2,1),Inf*ones(2,1));
phat
xSol
%% point #7 - displace, stretch, rotate
P = [0;0;0;1]; % <- this is a point in 3D affine space OR 4D
dx = [1 2 3];
DIS = [[1,0,0,dx(1)];[0,1,0,dx(2)];[0,0,1,dx(3)];[0,0,0,1]];
displacedP = DIS*P;
S = diag([1 2 3 1]);
DIS*S*P; % <- stretch and displace


%% use geometry for the fitting
% run PCA to get the expected parameters
[~,~,U,E,L,ERR,LAM] = PCA_FIT_FULL(d(:,1:2),2,1);
% get the angle
ang = atan2(E(2,1),E(1,1));
% make 2D hyperparameter distribution
[TM,iTM] = ndmap(2);
% fit mixture model with 1 distribution
GMModel = fitgmdist(d(:,1:2),1,'Options',ops);
% stack the data
data = [ones(size(d,1),1) d(:,1:2)]; 
% make the gaussian distribution
f = @(x,p)(prod(p(4:5))*(2*pi)^2)^-1*exp(-.5*sum(x(:,2:end).^2,2));
% make contrast function
func = @(p)-sum(log(f((iTM(p)*data')',p)));
% init parameter value
initP = [ang U (diag(LAM).^.5)'];
options = optimoptions('fmincon','Display','iter');
xSol = fmincon(func,initP,[],[],[],[],[],[]);
%xSol = fminsearch(func,initP);
% plot the reseults
ix = (iTM(xSol)*data')';close all;
plot(ix(:,2),ix(:,3),'.')
% the results check out
std(ix,1,1);
diag(LAM).^.5
xSol(4:5)
diag(GMModel.Sigma.^.5)
%% 
%% show how PCA is geometry
alpha = pi/5;
RM = [[cos(alpha) sin(alpha)];[-sin(alpha) cos(alpha)]];
RM = [[RM,[0;0]];[0 0 1]];
DIS = eye(3);
DIS(1:2,3) = [.3;1.4];
S = diag([3 1/2 1]);
T = DIS*RM*S;
[ty,tx] = ndgrid(linspace(-5,5,100),linspace(-5,5,100));
X = [tx(:),ty(:),ones(size(tx(:)))];
Y = (T*X')';
dx = makedist('norm',0,1);
dy = makedist('norm',0,1);
samx = random(dx,1000,1);
samy = random(dy,1000,1);
samTot = [samx samy ones(size(samx))];
samY = (T*samTot')'
f = mvnpdf(X(:,1:2),[0 0],eye(2));
f = reshape(f,size(tx));
close all
imshow(f,[]);
waitforbuttonpress
X = reshape(X,[size(tx) 3]);
Y = reshape(Y,[size(tx) 3]);
close all
figure;hold on
for e = 1:size(X,1)
    plot(X(e,:,1),X(e,:,2),'r');
end
for e = 1:size(X,2)
    plot(X(:,e,1),X(:,e,2),'r');
end
for e = 1:size(Y,1)
    plot(Y(e,:,1),Y(e,:,2),'b');
end
for e = 1:size(Y,2)
    plot(Y(:,e,1),Y(:,e,2),'b');
end
plot(samY(:,1),samY(:,2),'g.')
mean(samY)
DIS(1:2,3)
[S C U E L ERR LAM] = PCA_FIT_FULL(samY,2,1);
%% Test #2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FilePath = '/mnt/spaldingdata/nate/octerineDataStore/craineBrothers/fieldQuinoa/';
FileList = {};
FileExt = {'JPG'};
FileList = fdig(FilePath,FileList,FileExt,1);
%% gather stats on 5 files for 9 million points
Nfiles = 5;
Npoints = 10000;
[density] = channelSampler(FileList,Nfiles,Npoints);
%% fit n-gaussians
% try circle displace,rotated and streched
[TM,iTM] = ndmap(2);
%% make grid and transform
w = [pi/5 3 0 2 2];
[ty,tx] = ndgrid(linspace(-5,5,10),linspace(-5,5,10));
X = [ones(size(tx(:))),tx(:),ty(:)];
Y = (TM(w)*X')';
X = reshape(X,[size(tx) 3]);
Y = reshape(Y,[size(tx) 3]);
close all
figure;hold on
for e = 1:size(X,1)
    plot(X(e,:,2),X(e,:,3),'r');
end
for e = 1:size(X,2)
    plot(X(:,e,2),X(:,e,3),'r');
end
for e = 1:size(Y,1)
    plot(Y(e,:,2),Y(e,:,3),'b');
end
for e = 1:size(Y,2)
    plot(Y(:,e,2),Y(:,e,3),'b');
end
%% make master map

%%
szD = size(density.density);
nd = numel(szD);
P = randomPDF(3,[4 4 4]);
% make hypersphere
%%%%%%%%%%%%%%%%%%%%%
[f,TM,x] = hSphere(nd);
% make displacement
%%%%%%%%%%%%%%%%%%%%%
d = sym('d',[1 (nd)]);
dx = [];
for e = 1:nd
    dx = [dx ; symfun(d(e),[x d])];
end
TM = [[symfun(TM,[x d]),dx];[zeros(1,nd) 1]];
g = matlabFunction(f);

gTM = matlabFunction(inv(TM));
gTM = vectorizeFunc(gTM);

forwardMap = vectorizeFunc(matlabFunction(TM));
backwardMap = vectorizeFunc(matlabFunction(inv(TM)));
%%
x = [1 0 0 3 4 5];
forwardMap(x)
%%
gTM(p)
g = @(x,p)mtimesx(gTM(p),x);
p = [1 0 0 0 0 0];
x =[1;1;1];
%% prob map
close all
I = imread(FileList{2});
pI = density.prob(Io);
imshow(pI,[]);
%%
cos(p(1))
sin(p(1))*cos(p(2))
sin(p(1))*sin(p(2))
cos

pdf1 = @(x,p)mvnpdf(x,p(1:3),.5*(reshape(p(4:end),[3 3]) + reshape(p(4:end),[3 3])'));
.5*(reshape(p(4:end),[3 3]) + reshape(p(4:end),[3 3])')
pdf1([1 1 1],1:12)