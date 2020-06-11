%%
Q = 8;
phi1 = 2*pi*rand(Q,1) + pi;
phi2 = 2*pi*rand(Q,1) + pi;
mag0 = rand(Q,1);
para = [mag0 phi1 phi2];
reg = makeQuantumRegister(para(:));
prob = measureRegisterProb(reg);
close all
plot(prob)
%%
n = 4;
m = 4;
n = 256;
m = n;
xv = 0:(n-1);
yv = 0:(m-1);
m = yv'*xv;
M = (n^-.5)*exp((2*pi*1i*m)/n);

%% 

Q = 8;
phi1 = 2*pi*rand(Q,1) + pi;
phi2 = 2*pi*rand(Q,1) + pi;
mag0 = rand(Q,1);
para = [mag0 phi1 phi2];


close all
u = 6;
s = 0;
dm = (1:2^Q)';
pt = normpdf(dm,u,s);
clear x d
forceFunction = @(x)norm(pt - measureRegisterProb(makeQuantumRegister(x)));
forceFunction = @(x)norm((u - expectedValue(x,dm)));
expectedValue = @(x,d)d'*measureRegisterProb(makeQuantumRegister(x));
smoothness = @(x)max(abs(diff(x)));
d1 = @(x,d)norm(u - expectedValue(x,d));
d2 = @(x,d)expectedValue(x,(d-expectedValue(x,d)).^2);
d3 = @(x,d)quantumEntropy1(measureRegisterProb(makeQuantumRegister(x)));
mv = @(x,d)expectedValue(x,(d-expectedValue(x,d)).^2);
forceFunction = @(x)(-quantumEntropy1(measureRegisterProb(makeQuantumRegister(x))) + norm([u,s.^.5] - [expectedValue(x,dm) mv(x,dm).^.5]));
%forceFunction = @(x)(norm([u,s.^.5] - [expectedValue(x,dm) mv(x,dm).^.5]));
%forceFunction = @(x)(smoothness(x) + (norm([u,s.^.5] - [expectedValue(x,dm) mv(x,dm).^.5])));
%forceFunction = @(x)(d1(x,dm) + norm(pt - measureRegisterProb(makeQuantumRegister(x))));
%[e] = quantumEntropy1(p);
%forceFunction = @(x)norm([u] - [expectedValue(x,dm)]);
%forceFunction = @(x)norm([s] - [mv(x,dm)]);


forceFunction = @(x)(d1(x,dm)+(d2(x,dm).^.5)-d3(x,dm));
forceFunction = @(x)(d1(x,dm)-d3(x,dm));
forceFunction = @(x)(d1(x,dm)+(d2(x,dm).^.5));
contrast1 = @(x)forceFunction(x);


lb = [zeros(Q,1);-pi*ones(Q,1);-pi*ones(Q,1)];
ub = [ones(Q,1);pi*ones(Q,1);pi*ones(Q,1)];
xsol = fmincon(contrast1,para(:),[],[],[],[],lb,ub);
%xsol = patternsearch(contrast1,para(:),[],[],[],[],lb,ub);

p = measureRegisterProb(makeQuantumRegister(xsol));
plot(dm,p);hold all
qsol = M*makeQuantumRegister(xsol);
q = measureRegisterProb(qsol);
figure;plot(dm,q);
%%
contrast2 = @(x)norm([s] - [mv(x,dm)]);
noncon = @(x)[0 norm([u] - [expectedValue(x,dm)])];
QQQ = @(x)abc(noncon,x);
xsol = fmincon(contrast2,xsol,[],[],[],[],lb,ub,QQQ);
p = measureRegisterProb(makeQuantumRegister(xsol));
plot(p);
contrast1(xsol)
contrast2(xsol)
forceFunction = @(x)norm([u,s] - [expectedValue(x,dm) mv(x,dm)]);
forceFunction(xsol)

%%
qr = makeQuantumRegister(xsol)
%%
xSol = fminsearch(contrast,para(:));
reshape(xSol,[8 3])
