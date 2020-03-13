C = readtext('/home/nate/ABCD for MLM.csv');
C(:,1) = [];
%{
for e = 1:size(C,2)
    C{1,e} = strrep(C{1,e},'Col-0','col0');
end
%}
%%
D = cell2mat(C(3:end,1:end));
[cS cC cU cE cL cERR cLAM] = PCA_FIT_FULL_T(D,1);
tbl = table(C(1,1:end)',cell2mat(C(2,1:end)'),cC(1,:)','VariableNames',{'Genotype','Stress','PC1'});
%%
t = table;
for e = 1:size(C,2)
    NN{e} = [C{1,e} '-stress' num2str(C{2,e})];
end
UQ = unique(NN');
eq = ['PC1~-1+'];
int = '';
EE = 46;
for u = 1:numel(UQ)
    IDXE{u} = strcmp(NN',UQ{u});
    UQ{u} = strrep(strrep(lower(UQ{u}),'_',''),'-','');
    t{:,u} = (IDXE{u});
    if u < EE
        eq = [eq UQ{u} '+' ];
        int = [int '+' UQ{u} ':stress'];
    end
    
end
eq(end) = [];
%t{:,end+1} = cell2mat(C(2,1:end)');
t{:,end+1} = cC(1,:)';
%UQ{end+1} = 'stress';
UQ{end+1} = 'PC1';
t.Properties.VariableNames = UQ;
lme = fitlme(t,[eq]);
%%
lme = fitlme(tbl,'PC1~1+Genotype+Stress+Genotype:Stress');
%%
fidx = find(strcmp(C(1,:),'Col-0') & (cell2mat(C(2,:)) == 0));
fidx1 = find(strcmp(C(1,:),'cam2_3_5') & (cell2mat(C(2,:)) == 1));
%fidx = find(strcmp(C(1,:),'cam7') & (cell2mat(C(2,:)) == 0));
levelM = mean(cC(:,fidx))
%%cC = bsxfun(@minus,cC,mean(cC(:,fidx)));
[h,p] = ttest2(cC(fidx),cC(fidx1));
%%
%lme = fitlme(tbl,'PC1~1+(Genotype|Stress)');
%lme = fitlme(tbl,'PC1~1+(Genotype|Stress)');
lme = fitlme(tbl,'PC1~1+Genotype+Stress+Genotype:Stress');
G = lme.predict;
mean(G(fidx))

%%
[B,Bnames,stats] = randomEffects(lme)
%%
base = [0 zeros(1,21),0,zeros(1,21)];
%base = 0;
n = base;
base(14) = 1;
%base(23) = 1;
%base(1) = 1;
%base(5) = 1;
%base(6) = 1;
base(18) = -1;
%base(40) = -1;

%
H = [base;n];
%K = zeros(1,length(stats));
[pVal,F,DF1,DF2] = coefTest(lme,base)




%%
fidxMT = find(strcmp(stats.Level,'cam2_3_5'));
fidxWT = find(strcmp(stats.Level,'col0'));
stress = 2;
nstress = 1;
K = zeros(1,length(stats));
K(1,fidxWT(stress)) = 1;
coefTest(lme,base,levelM,'REContrast',K)
%%
close all
G = lme.predict;
SIM = PCA_BKPROJ_T(G',cE,cU);

for tr = 1:size(SIM,2)
    plot(cS(:,tr),'b');
    hold on
    plot(SIM(:,tr),'r')
    hold on
    plot(D(:,tr),'k')
    hold off
    %waitforbuttonpress
end
%%
toInclude = zeros(size(lme.PredictorNames));
toInclude(42) =  1;
toInclude(40) = 1;
tblN = table;
for e =1:44
    tblN(1,lme.PredictorNames{e}) = {logical(toInclude(e))};
end

GG = PCA_BKPROJ_T(lme.predict(tblN),cE,cU);
%%
close all
GENO = unique(C(1,2:end));
for u = 1:numel(GENO)
    close all
    
    
    
    %GENO{u} = 'cml49'
    fidx = find(strcmp(C(1,2:end),GENO{u}))
    stress = cell2mat(C(2,2:end));
    g0 = find(stress(fidx)==0);
    g1 = find(stress(fidx)==1);
    sub = C(4:end,fidx(g0));
    U0 = mean(cell2mat(sub),2);
    S0 = std(cell2mat(sub),1,2)*size(sub,2)^-.5;
    sub = C(4:end,fidx(g1));
    U2 = mean(cell2mat(sub),2);
    S2 = std(cell2mat(sub),1,2)*size(sub,2)^-.5;
    
    
    
    
    
    
    
    
    errorbar(U0,S0,'b');
    
    hold all
    errorbar(U2,S2,'c');
    
    
    
    
    
    
    hold all
    
    
    
    
    
    
    %{
    
    
    SSS = std(SIM(:,fidx(g0)),1,2)*numel(g0)^-.5;
    UUU = mean(SIM(:,fidx(g0)),2);
    errorbar(UUU,SSS,'k');
    
    
    SSS = std(SIM(:,fidx(g1)),1,2)*numel(g1)^-.5;
    UUU = mean(SIM(:,fidx(g1)),2);
    errorbar(UUU,SSS,'c');
    hold on
    title(GENO{u})
    
    %}
    
    
    
    
    
    fidx = find(strcmp(C(1,2:end),'Col-0'))
    stress = cell2mat(C(2,2:end));
    g0 = find(stress(fidx)==0);
    g1 = find(stress(fidx)==1);
    sub = C(4:end,fidx(g0));
    U0 = mean(cell2mat(sub),2);
    S0 = std(cell2mat(sub),1,2)*size(sub,2)^-.5;
    sub = C(4:end,fidx(g1));
    U2 = mean(cell2mat(sub),2);
    S2 = std(cell2mat(sub),1,2)*size(sub,2)^-.5;
    
    
    
    
    errorbar(U0,S0,'k');
    
    hold all
    errorbar(U2,S2,'k'); 
    title(GENO{u});
    
    plot(GG,'r')
    waitforbuttonpress
end