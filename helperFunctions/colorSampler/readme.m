% readme for list of helper function(s)
% channelSampler - sample color(or other) channel(s) to build up
% distribution. The color distribution will be a dataTensor.
%% date: April 13, 2020
% attempting to use this with sampling the velocity and intensity in
% shistosome datasets.  The reader is now a handle which can be passed in.
% The reader handle takes the file to read and a number which is the number
% of the file in a list.  This is used to load assocated data, like a mask,
% from a list.
%% date: April 16, 2020
% Julian and I are testing the discrete distribution via pc scores from
% image patches
% data @ '/mnt/spaldingdata/nate/octerineDataStore/testData/channelSampler/200416_SegmentScores_3879Patches.mat'
%% date: April 28, 2020
% Joe and I needed the project to work. AKA: it would have made it easier
% for programming a solution.
%%  
%%
%project.save('colorSampler')
project.save()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% quick test
r1 = random('norm',40,100,1000,1);
r2 = random('norm',20,3,1000,1);
r3 = random('norm',9,6,1000,1);
data = [r1 r2 r3];
test_pdf = @(d,x)mvnpdf(d,x(1:3),diag(x(4:6)));
test_pdf(data,rand(1,6))
%% this does not work because mle is vector only
% therefore the below search is needed
phat = mle(data,'pdf',test_pdf,'start',rand(1,6));
%%
N = 10000;
r1 = random('norm',40,100,N,1);
r2 = random('norm',20,3,N,1);
r3 = random('norm',9,6,N,1);
data = [r1 r2 r3];
test_pdf = @(d,x)mvnpdf(d,x(1:3),abs(diag(x(4:6))));
mag = 5;
p = [40 20 9 100 3 6];
pr = p + mag*rand(1,6);
myLogContrast = @(x,p,func)-sum(log(func(x,p)));
myLogContrast(data,pr,test_pdf);
clear para


func = @(para)myLogContrast(data,para,test_pdf);
contrast = func;
initX = pr;

solverType = 'patternsearch';
%solverType = 'fminsearch';
%solverType = 's';
%jumpMAX = 15;
switch solverType
    case 'fminsearch'
        options = optimset('Display','iter','TolX',10^-6,'TolFun',10^-6);
        xSol = fminsearch(contrast,initX,options);
    case 'patternsearch'
        A = [];b = [];Aeq = [];beq = [];lb = -ones(numel(initX),1);ub = -lb;nonlcon = [];
        A = [];b = [];Aeq = [];beq = [];lb = zeros(numel(initX),1);ub = ones(numel(initX),1);nonlcon = [];
        %A = [];b = [];Aeq = [];beq = [];lb = [];ub = [];nonlcon = [];
        options = optimoptions('patternsearch','Display','none');
        xSol = patternsearch(contrast,initX,A,b,Aeq,beq,lb,ub,nonlcon,options);
    case 'fmincon'
        A = [];b = [];Aeq = [];beq = [];lb = -ones(numel(initX),1);ub = -lb;nonlcon = [];
        options = optimoptions('fmincon','Display','iter','FiniteDifferenceType','central',...
            'StepTolerance',10^-10);
        xSol = fmincon(contrast,initX,A,b,Aeq,beq,lb,ub,nonlcon,options);
    case 's'
        options = optimset('Display','iter','TolX',10^-6,'TolFun',10^-6);
        [xSol1,v1] = fminsearch(contrast,initX,options);
        options = optimoptions('patternsearch','Display','iter','TolMesh',10^-10);
        %A = [];b = [];Aeq = [];beq = [];lb = zeros(numel(initX),1);ub = ones(numel(initX),1);nonlcon = [];
        A = [];b = [];Aeq = [];beq = [];lb = [-jumpMAX -jumpMAX -pi/2 .9 .9 -pi/8];ub = [jumpMAX jumpMAX pi/2 1.1 1.1 pi/8];nonlcon = [];

        [xSol,v2] = patternsearch(contrast,xSol1,A,b,Aeq,beq,lb,ub,nonlcon,options);

        v1
        v2
        %{
        if v2 <= v1
            xSol = xSol2;
        else
            xSol = xSol1;
        end
        %}
end