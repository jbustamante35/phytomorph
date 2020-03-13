function model=PCAGCA(X,options)
% options=PCAGCA(X)
% model=PCAGCA(X,options)
%
% Separates common and distinctive components in multiple data blocks,
% using a combination of PCA and GCA. reference: Smilde AK, M�ge I, N�s T, et al. Common and Distinct Components in Data Fusion. J Chemom. 2017;In press. doi:10.1002/cem.2900.
%
%
% X is cell array of "saisir" data sets. The "saisir" data set structure
% contains three fields:
%     .d: The data matrix, size NxP
%     .i  Sample names (character array, N rows)
%     .v  Variable names (character array, P rows)
% (See http://www.chimiometrie.fr/saisir_webpage.html for more information about the saisir toolbox)
%
%
% options=PCAGCA(X) gives default options, which is a struct with fields:
%  -BlockNames: character array with rows euql to number of blocks
%  -commondim: either "rows" or "columns"
%  -nObj: Number of objects (in each block)
%  -nVars: Number of variables (in each block)
%  -autoselect: 0/1 (interactive or automatic selection of components.
%           Interactive is recommended.
%  -Rlim: Threshold for correlation when common components are selected
%           automatically (only used if if autoselect=1)
% -ExpVarLim: Threshold for explained variance (per block) when common components are selected
%           automatically (only used if autoselect=1)
% -preproc: cell array (length nBlock) with preprocessing option for each
%           block. Each cell should be either 'mean center','autoscale' or
%           'none'.
% -Amax:    Vector (size nBlocks x 1). Maximum number of components for
%           each block.
% -blockCombinations: cell array, the cells define in which block
%           combinations to search for common components, as well as the order. Each
%           cell is a vector with block indices.
% -nCompsForEachContribution: cell array, same size as "blockCombinations".
%           Defines the number of common components for each subspace.
% -nCompsLocalPCA: vector(size nBlocks x 1). Defines the number of
%           components to extract from each individual PCA model, to be
%           used for the canonical correlation analysis.
%
%
%
% model=PCAGCA(X,options) calculates model parameters. "model" is a struct
% with fields:
%     -X: The input data (cell array)
%     -options: same as the options struct described above
%     -Aselected: (vector, length nBlocks). Contains the number of components that have been extracted from each individaul PCA model
%     -Scores: cell array (length nBlcoks) containing the scores for each block
%     -Loadings: cell array (length nBlcoks) containing the loadings for each block
%     -ExpVar: cell array (length nBlcoks) containing the component-wise explained variances for each block
%     -R: vector (length equal to the total number of common components) containing the canonical correlations.
%     -Labels: cell array (length nBlcoks) containing the component labels for each block. Prefix 'C' stands for Common and 'D' stands for distinct
%
% Written by Ingrid M�ge, modified August 2017.

%check if rows or columns are common dimension


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==1
    
    nBlocks=length(X);
    model.BlockNames=strcat('Block-',num2str((1:nBlocks)'));
    
    %chek common dimension (rows or columns)
    for i=1:nBlocks
    [n(i), d(i)] =size(X{i}.d);
    end
    if std(n)==0,
    model.commondim='rows'; % rowwise linkage.. same objects
    model.nObj=n(1);
    model.nVars=d;
    else
    model.commondim='columns'; % variablewise linkage.. same objects  
    model.nObj=n;
    model.nVars=d(1);
    end
    
    model.autoselect=0;
    model.Rlim=[];
    model.ExpVarLim=[];
    for i=1:nBlocks;
        model.preproc{i}='mean center';
        model.Amax(i)=min(rank(X{i}.d)-1,20);
    end
    
    
    model.blockCombinations={};
    for i=1:(length(X)-1)
        r=combnk(1:length(X),length(X)-i+1);
        n=size(r,1);
        model.blockCombinations(end+1:end+n)=mat2cell(r,ones(n,1),length(X)-i+1)';
    end
    model.nCompsForEachContribution=cell(1,length(model.blockCombinations));
    model.nCompsLocalPCA=[];
    model.BlockNames=strcat('Block-',num2str((1:nBlocks)'));
    
    
    
else
    
    
    nBlocks=length(X);
    nObj=options.nObj;
    nVars=options.nVars;
    autoselect=options.autoselect;
    Amax=options.Amax;
    blockCombinations=options.blockCombinations;
    nCompsForEachContr=options.nCompsForEachContribution;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initialise model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    model.type = 'PCA-GCA';
    model.X=X; %save raw data
    model.options=options;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % preprocess X and Y
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i=1:nBlocks
        X{i}=mypreprocess(X{i},options.preproc{i});
        X{i}=rmfield(X{i},'pp');
        totVar{i}=sum(sum(X{i}.d.^2));
        Labels{i}=char();
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initial PCA for each block
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i=1:nBlocks
        
        m(i)=mypca(X{i},Amax(i),options.preproc{i});
        if isempty(options.nCompsLocalPCA)
        if autoselect==1
            A(i)=m(i).Aopt;
        else
            figure; pcaplots(m(i));
            A(i)=input(['Number of components in block ' num2str(i) ': ']);
        end
        else
           A(i)=options.nCompsLocalPCA(i);
        end
        
        T{i}=mat2saisir(m(i).T(:,1:A(i)));
        P{i}=mat2saisir(m(i).P(:,1:A(i)));
        Scores{i}=[];
        ExpVar{i}=[];
        Loadings{i}=[];
        
    end
    model.Aselected=A;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform CCA on each of the blockCombinations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Scores=cell(1,length(options.blockCombinations));
    %Loadings=cell(1,length(options.blockCombinations));
    
    R=[];
    
    for k=1:length(blockCombinations)
        idx=blockCombinations{k}; %which blocks to search for common components in
        
        
        
        %Canonical correlation analysis
        if strcmp(options.commondim,'rows')==1 %rows are common dimension, GCA is run on scores T
        [A,r,U,TC]=GCCA(T(idx));
        else %columns are common dimension, GCA is run on the loadings P
        [A,r,U,TC]=GCCA(P(idx));
        end
   
        
        if isempty(options.nCompsForEachContribution{k}) %the number of components is not pre-defined
        %calculate explained variances for each canonical component
        expVarX = zeros(length(r),length(idx));
        for ii=1:length(idx)
            XX=X{idx(ii)}.d;
            for jj=1:length(r)
                ss=U{ii}(:,jj)/sqrt(U{ii}(:,jj)'*U{ii}(:,jj));
                if strcmp(options.commondim,'rows')==1 %rows are common dimension
                loads=XX'*ss;
                XX=XX-ss*loads';
                else
                scores=XX*ss;
                XX=XX-scores*ss';
                end
                expVarX(jj,ii)=(1-sum(sum(XX.^2))/totVar{idx(ii)})*100;
            end
        end
        expVarX=diff([repmat(0,1,length(idx));expVarX]);
        % plot explained variances and correlations, select components
        
        if options.autoselect==0
        h=figure;
        set(h,'Name',['Common components, block ' num2str(idx)])
        if length(r)==1
            hh=bar([expVarX; repmat(0,1,length(idx))]);
        else
            hh=bar([expVarX]);
        end
        hold on
        
            legend(options.BlockNames(idx,:));
     
        
        plot(1:length(r),r*100,'.','markersize',20,'color','b')
        line([0 length(r)+1],[90 90],'linestyle',':')
        line([0 length(r)+1],[95 95],'linestyle',':')
        title('R and explained variances for each canonical component')
        %save numbers for plot
        model.numbers_for_plot{k}=[r' expVarX];
        
        
        nComps = input(['Common components, block ' num2str(idx) ', keep which components (vector)? ']);
        if length(nComps)==1; nComps=1:nComps; end
        
        else %automatic component selection based on pre-defined limits on R and ExpVar 
            nComps =find(r>options.Rlim & min(expVarX')>options.ExpVarLim);
            model.options.nCompsForEachContribution{k}=nComps;
        end
            
        else nComps=options.nCompsForEachContribution{k}; %the number of components were set in options
        end
        
        if nComps==0; nComps=[]; end
        model.options.nCompsForEachContribution{k}=nComps;
        %calculate scores/loadings and orthogonalize T/P
        if ~isempty(nComps)
            R(end+1:end+length(nComps))=r(nComps);
            
            for iblock=1:length(idx)
                if strcmp(options.commondim,'rows')==1 %rows are common dimension
                Scores{idx(iblock)}=[Scores{idx(iblock)} U{iblock}(:,nComps)./repmat(sqrt(diag(U{iblock}(:,nComps)'*U{iblock}(:,nComps)))',nObj,1)];
                if nBlocks==2; Labels{idx(iblock)}=strcat('C_',num2str((1:length(nComps))')); else Labels{idx(iblock)}=strvcat(Labels{idx(iblock)},strcat(strcat('C', num2str(idx')','_'),num2str((1:length(nComps))'))); end
                T{idx(iblock)}.d=T{idx(iblock)}.d-U{iblock}(:,nComps)*inv(U{iblock}(:,nComps)'*U{iblock}(:,nComps))*U{iblock}(:,nComps)'*T{idx(iblock)}.d;
                
                else %columns are common dimension
                                  Loadings{idx(iblock)}=[Loadings{idx(iblock)} U{iblock}(:,nComps)./repmat(sqrt(diag(U{iblock}(:,nComps)'*U{iblock}(:,nComps)))',nVars,1)];
                if nBlocks==2; Labels{idx(iblock)}=strcat('C_',num2str((1:length(nComps))')); else Labels{idx(iblock)}=strvcat(Labels{idx(iblock)},strcat(strcat('C', num2str(idx')','_'),num2str((1:length(nComps))'))); end
                P{idx(iblock)}.d=P{idx(iblock)}.d-U{iblock}(:,nComps)*inv(U{iblock}(:,nComps)'*U{iblock}(:,nComps))*U{iblock}(:,nComps)'*P{idx(iblock)}.d;
                end

            end
        end
    end
    
    %Restructure unique part in each block
    
    for iblock=1:nBlocks
        
        
        if strcmp(options.commondim,'rows')==1 %rows are common dimension
        %check if more dimension in X
        dim_left=model.Aselected(iblock)-size(Scores{iblock},2);
        if dim_left>0
            [uu,ss,vv]=svds(T{iblock}.d,dim_left);
            Scores{iblock}=[Scores{iblock} uu];
            Labels{iblock}=strvcat(Labels{iblock},strcat('D',num2str(iblock),'_',num2str((1:dim_left)')));
        end
        
        else %columns common dimension
        %check if more dimension in X
        dim_left=model.Aselected(iblock)-size(Loadings{iblock},2);
        if dim_left>0
            [uu,ss,vv]=svds(P{iblock}.d,dim_left);
            Loadings{iblock}=[Loadings{iblock} uu];
            Labels{iblock}=strvcat(Labels{iblock},strcat('D',num2str(iblock),'_',num2str((1:dim_left)')));
        end

            
        end
      
        
    end
    
    %calculate scores/loadings for all blocks and components
    for iblock=1:nBlocks
        if strcmp(options.commondim,'rows')==1 %rows are common dimension
        for icomp=1:size(Scores{iblock},2)
            Loadings{iblock}(:,icomp)=X{iblock}.d'*Scores{iblock}(:,icomp);
            X{iblock}.d=X{iblock}.d-Scores{iblock}(:,icomp)*Loadings{iblock}(:,icomp)';
            ExpVar{iblock}(icomp)=(1-sum(sum(X{iblock}.d.^2))/totVar{iblock})*100;
        end
        
        else %columns common dimension
             for icomp=1:size(Loadings{iblock},2)
            Scores{iblock}(:,icomp)=X{iblock}.d*Loadings{iblock}(:,icomp);
            X{iblock}.d=X{iblock}.d-Scores{iblock}(:,icomp)*Loadings{iblock}(:,icomp)';
            ExpVar{iblock}(icomp)=(1-sum(sum(X{iblock}.d.^2))/totVar{iblock})*100;
             end
        end
            
         %remove small distinct components if threshold for ExpVar is set
        if ~isempty(options.ExpVarLim)
            delidx=find(diff([0 ExpVar{iblock}])<options.ExpVarLim);
            Scores{iblock}(:,delidx)=[];
            Loadings{iblock}(:,delidx)=[];
            ExpVar{iblock}(delidx)=[];
            Labels{iblock}(delidx,:)=[];
        end
    end
    
     
    
    model.Scores=Scores;
    model.Loadings=Loadings;
    for i=1:nBlocks; model.ExpVar{i}=diff([0 ExpVar{i}]); end
    model.R=R;
    model.Labels=Labels;
    
    model.ResTab=results_table(model);
    
end

function X=mypreprocess(X,method)
% X=mypreprocess(X,method)
%
% method is a string, either 'mean center' or 'autoscale'

if isempty(method)
    X.pp.meth='No preprocessing';
    X.pp.mean=zeros(1,size(X.d,2));
    X.pp.std=ones(1,size(X.d,2));
else


        if isequal(lower(method),'auto scale')
            method='autoscale';
        end
        switch lower(method)

            case 'mean center'
                if isfield(X,'pp')
                    if isequal(lower(X.pp.meth),'autoscale') | isequal(lower(X.pp.meth),'mean center')
                        error('Data set is alreay centered!')
                    end
                end
                Xm=mean(X.d);
                Xs=ones(1,size(X.d,2));
                X.d=(X.d-repmat(Xm,size(X.d,1),1))./repmat(Xs,size(X.d,1),1);
                X.pp.meth=method;
                X.pp.mean=Xm;
                X.pp.std=Xs;


            case 'autoscale'
                if isfield(X,'pp')
                    if isequal(lower(X.pp.meth),'autoscale')
                        error('Data set is alreay autoscaled!')
                    end
                end
                Xm=mean(X.d);
                Xs=std(X.d);
                X.d=(X.d-repmat(Xm,size(X.d,1),1))./repmat(Xs,size(X.d,1),1);
                X.pp.meth=method;
                X.pp.mean=Xm;
                X.pp.std=Xs;
            otherwise
                if ~isfield(X,'pp')
                    X.pp.meth='No preprocessing';
                    X.pp.mean=zeros(1,size(X.d,2));
                    X.pp.std=ones(1,size(X.d,2));
                end
        end
end

function model=mypca(X,A,preproc)
% model=mypca(X,A,preproc)
% 
% 
% INPUT
% X         Descriptor data matrix, saisir structure
% A         Maximum number of components
% preproc  optional. cell array containing preprocessing methods for X
% 
% OUTPUT
% model
% 
% Created: 10/11 2008
% Modified: 
% Status:  works 10/11 2008
%
% Ingrid M�ge


[n,p]=size(X.d);
nx = n;
px = p;

%set labels if they are not given in saisir
if  ~isfield(X,'v')
    X.v=num2str((1:px)');
elseif isempty(X.v)
    X.v=num2str((1:px)');
end

if ~isfield(X,'i')
    X.i=num2str((1:nx)');
elseif isempty(X.i) 
    X.i=num2str((1:nx)');
end

model.type = 'pca';
model.X=X; %save raw data

if nargin==3
    X=mypreprocess(X,preproc);
    model.Xpp=X.pp;
    X=rmfield(X,'pp');
end


[U,S,P]=svd(X.d);
T=U*S;
T=T(:,1:A);
P=P(:,1:A);
totvar=sum((diag(S).^2));
S=S(1:A,1:A);

ExpVar=[0; cumsum(diag(S.^2)./repmat(totvar,A,1)*100)]';

E=cell(A+1,1);
ResVar=zeros(A+1,1);
ResVarVal=zeros(A+1,1);
SampRes=zeros(n,A+1);
VarRes=zeros(p,A+1);
Hi=zeros(n,A+1);
for i=0:A
    E{i+1}=X.d-T(:,1:i)*P(:,1:i)';
    ResVar(i+1)=mean(mean(E{i+1}.^2));
    SampRes(:,i+1)=mean((E{i+1}.^2)')';
    VarRes(:,i+1)=mean((E{i+1}.^2))';
    if i>0
        Hi(:,i+1)=Hi(:,i)+T(:,i).^2/(T(:,i)'*T(:,i));
    end
    ResVarVal(i+1)=mean(mean((E{i+1}).^2./repmat((1-Hi(:,i+1)).^2,1,p)));

end
   
Hi(:,1)=[];
Hi=Hi+1/n;

model.T=T;
model.P=P;
model.ResVar=ResVar;
model.ResVarVal=ResVarVal;
model.ExpVar=ExpVar;
model.SampRes=SampRes;
model.VarRes=VarRes;
model.VarExp=100-VarRes'./repmat(VarRes(:,1)',A+1,1)*100;
model.Hi=Hi;

crit = (1:A)'*0.02*mean(ResVar(1)) + mean(ResVar(2:end),2);
[dummy, Aopt] = min(crit);
model.Aopt = Aopt;
model.LevLim=3;
model.p=0.95;


function []=pcaplots(model,A1,A2,plottype,spec)
% []=pcaplots(model,A1,A2,plottype)
%
%
% Plots:
% scores for A1 and A2,
% loadings for A1 and A2
% explained variance for all components
% normal probability plot of residuals
%
% INPUT
% model       model struct from mypls
% A1,A2       PCs to plot for scores and loadings
% plottype    Optional.   =1: plots one figure with subplots (default)
%                         =2: plots each plot in a separate figure
% spec        Optional.   =1: plots loadings as curves

if ~exist('A1','var')
    A1=1;
end

if ~exist('A2','var')
    A2=2;
    if size(model.T,2)==1; A2=1; end
end

if ~exist('plottype','var')
    plottype=1;
end

if ~exist('spec','var')
    spec=0;
end



%Scores
if plottype==1
    figure('Name','PCA overview');
subplot(2,2,1)
else
    figure('Name','Scores');
end
H=plot(model.T(:,A1),model.T(:,A2),'.');
hold on
T1range=max(model.T(:,A1))-min(model.T(:,A1));
text(model.T(:,A1)+T1range/100,model.T(:,A2),model.X.i)
title('Scores')
xlabel(['PC ' num2str(A1) ', ' num2str(round(model.ExpVar(A1+1)-model.ExpVar(A1))) '%'])
ylabel(['PC ' num2str(A2) ', ' num2str(round(model.ExpVar(A2+1)-model.ExpVar(A2))) '%'])
grid on
h=gca;
if get(h,'Xlim')==0
    set(h,'Xtick',0)
else
set(h,'Xtick',sort(unique([get(h,'Xlim') 0])))
end
if get(h,'Ylim')==0
    set(h,'Ytick',0)
else
set(h,'Ytick',sort([get(h,'Ylim') 0]))
end


%Correlation Loadings
if plottype==1
    subplot(2,2,2)
else
    figure('Name','Loadings');
end

if spec~=1
hold off
plot(corr(model.T(:,A1),model.X.d),corr(model.T(:,A2),model.X.d),'.')
hold on
H=text(corr(model.T(:,A1),model.X.d)+2/100,corr(model.T(:,A2),model.X.d),model.X.v);
set(H,'Color','b')
title('Correlation Loadings')
xlabel(['PC ' num2str(A1)])
ylabel(['PC ' num2str(A2)])
%circle 100%
[X,Y] = pol2cart(linspace(0,2*pi,100),ones(1,100));
plot(X,Y,'k');
%circle 50%
[X,Y] = pol2cart(linspace(0,2*pi,100),ones(1,100)*sqrt(0.5));
plot(X,Y,'k');
H=line([0 0],[-1 1]);
set(H,'linestyle',':','Color','k')
H=line([-1 1],[0 0]);
set(H,'linestyle',':','Color','k')
else
    plot(str2num(model.X.v),model.P(:,1:model.Aopt)); axis tight
    legend(strcat('PC',num2str((1:model.Aopt)')))
    title('Loadings')
end

%Explained variance
if plottype==1
    subplot(2,2,3)
else
    figure('Name','Explained variance');
end


plot(0:size(model.T,2),model.ExpVar)
hold on
plot(model.Aopt,model.ExpVar(model.Aopt+1),'o')
title('Explained variance, Total')
xlabel('Components')
H=line([0 size(model.T,2)],[0 0]);
set(H,'color','k')

% leverage plot
if plottype==1
    subplot(2,2,4)
else
    figure('Name','Leverage');
end

plot(model.Hi(:,model.Aopt),model.SampRes(:,model.Aopt+1),'.')
Hirange=max(model.Hi(:,model.Aopt))-min(model.Hi(:,model.Aopt));
text(model.Hi(:,model.Aopt)+Hirange/100,model.SampRes(:,model.Aopt+1),model.X.i)
title(['Influence at ' num2str(model.Aopt) ' PCs'])
xlabel('Leverage')
ylabel('Residual')


function [A,R,U,C,predparam]=GCCA(X)
% [A,R,U,C,predparam]=GCCA(X)
% generalized canonical correlation analysis
% 
% INPUT
% X         Cell array of data saisir matrices
% 
% OUTPUT
% A         Cell array of canonical coefficients
% R         average corr between all pairs of canonical scores
% U         Cell array of canonical scores
% C         Matrix of consensus scores
% predparam parameters needed for prediction

ntable=length(X);

T={};
minrank=min(size(X{1}.d));
for i=1:ntable
    predparam.Mean{i}=mean(X{i}.d);
    X{i}.d=X{i}.d-repmat(mean(X{i}.d),size(X{i}.d,1),1);
    [U,S,V]=svd(X{i}.d);
    if size(S,2)>1
        thisrank=rank(S);
    else
        thisrank=1;
    end
    
    if thisrank<minrank
        minrank=thisrank;
    end
    predparam.Std{i}=std(U(:,1:thisrank));
    T{i}=U(:,1:thisrank)./repmat(std(U(:,1:thisrank)),size(X{i}.d,1),1);
    
    predparam.Pblock{i}=S(1:thisrank,1:thisrank)*V(:,1:thisrank)';
end
    
[C,S,V]=svd(cell2mat(T));
C=C(:,1:minrank);
predparam.Ptotal=S(1:minrank,1:minrank)*V(:,1:minrank)';

U={};
A={};
for i=1:ntable
    A{i}=pinv(X{i}.d'*X{i}.d)*X{i}.d'*C;
    U{i}=X{i}.d*A{i};
end

% lambda=zeros(1,minrank);
% for j=1:minrank
% for i=1:ntable
%     lambda(j)=lambda(j)+(corr(C(:,j),U{i}(:,j)))^2;
% end
% end
% lambda=lambda/ntable

ii=0;
R=zeros(1,minrank);
for i=1:ntable
    for j=i+1:ntable
ii=ii+1;
R=R+diag(corr(U{i},U{j}))';
    
end
end
R=R/ii;
  
function XX=mat2saisir(X,ObjLabels,VarLabels)

XX.d=X;

if exist('ObjLabels','var')
    XX.i=ObjLabels;
else
    XX.i=num2str((1:size(X,1))');
end

if exist('VarLabels','var')
    XX.v=VarLabels;
else
    XX.v=num2str((1:size(X,2))');
end

function ResTab=results_table(model)


allcomps=unique(char(model.Labels'),'rows');
nComps=size(allcomps,1);
nBlocks=length(model.Labels);

ResTab=cell(nComps+1,nBlocks+1);
for i=1:nComps
    for j=1:nBlocks
    
        idx=find(strcmp(strtrim(allcomps(i,:)),cellstr(model.Labels{j})));
        ResTab{i,j}=num2str(round(model.ExpVar{j}(idx)*10)/10); 
        ResTab{end,j}=num2str(round(sum(model.ExpVar{j})*10)/10);
    end
    
    if allcomps(i,1)=='C'; ResTab{i,nBlocks+1}=num2str(round(model.R(i)*100)/100); end
    
end

AttrNames=cell(1,nBlocks+1);
AttrNames(1:nBlocks)=cellstr(model.options.BlockNames)';
AttrNames{end}='Correlation';
%remove signs from attrnames
signs={'-','.',' ',':'};
for i=1:nBlocks
    for j=1:length(signs)
   idx=strfind(AttrNames{i},signs{j});
   AttrNames{i}(idx)=[];
    end
    if ~isvarname(AttrNames{i}); AttrNames{i}=['Block' num2str(i)]; end
end
Rownames=cellstr(allcomps); Rownames(end+1)={'Total'};
ResTab=cell2table(ResTab,'VariableNames',AttrNames,'RowNames',Rownames);

figure;
uitable('Data',ResTab{:,:},'ColumnName',ResTab.Properties.VariableNames,...
    'RowName',ResTab.Properties.RowNames,'Units', 'Normalized');
title('Explained variance per data block'); axis off



