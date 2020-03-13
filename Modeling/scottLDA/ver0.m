D = readtext('/mnt/spaldingdata/nate/GWAS_panel_minDP-10_max-missing-0.5-notimputed.hmp.txt');
fprintf(['Done loading data.\n'])
%% scan the data line by line to creat the genotype data
for line = 1:size(D,1)
   
    col = textscan(D{line},'%s','Delimiter','\t');
    col = col{1}';
    if line == 1
        T = cell(size(D,1),size(col,2));
    end
    for c = 1:numel(col)
        T{line,c} = col{c};
    end
    line
    size(D,1)
end
fprintf(['Done scanning data.\n'])
%%
tic
T = cellfun(@(x) myTextScan1(x), D, 'UniformOutput', 0);
toc
fprintf(['Done scanning lines.\n'])
%% remove the one row that does not match in szei
for e = 1:size(T,1)
    if e == 1
        sz = size(T{e});
    end
    ch(e) = all(size(T{e}) == sz);
end
ridx = find(ch==0);
T(ridx) = [];
% cat the cell array
tic
T = cat(1, T{:});
toc
fprintf(['Done stacking lines.\n'])
%% read the pheotype data
P = readtext('/mnt/spaldingdata/nate/phenotypes.csv');
%% read the key file
keyData = readtext('/home/nate/keyfile.csv');
%% chrom number
for e = 2:size(T,1)
    fidx = strfind(T{e,3},'_');
    chr(e) = str2num(T{e,3}((fidx(1)+4):end));
end
%% extract phenotype data
M = [];
uM = [];
s1 = {'A','T','G','C'};
s2 = {'A','T','G','C'};
[snpV] = scartprod(s1,s2);
for k = 1:size(keyData,1)
    tmpK1 = keyData{k,2};
    tmpK2 = keyData{k,4};
    tmpK1;
    tmpK2;
    class(tmpK1);
    
    tidx1 = find(strcmp(T(1,:),num2str(tmpK1)));
    tidx2 = find(strcmp(P(:,5),['"' tmpK2 '"']));
    
    % for the whole
    %tpM = P(tidx2,20:24);
    % for shoulder
    tpM = P(tidx2,25:29);
    
    
    if ~isempty(tpM)
        
        
        
        for e = 1:numel(tpM)
            tmp = tpM{e};
            tmp(1) = [];
            tmp(end) = [];
            tmp = str2num(tmp);
            tpM{e} = tmp;
        end
        tpM = cell2mat(tpM);
        
        
        
        tgM = T(2:end,tidx1);




        parfor s = 1:numel(tgM)
            gidx = find(strcmp(tgM{s},snpV));
            if isempty(gidx)
                tgM{s} = NaN;
            else
                tgM{s} = gidx;
            end
        end
        tgM = cell2mat(tgM)';
        
        
         
        UtpM = mean(tpM,1);
        UvecM = [k*ones(size(UtpM,1),1) UtpM tgM];
        
       
        tgM = repmat(tgM,[size(tpM,1) 1]);
        vecM = [k*ones(size(tpM,1),1) tpM tgM];
       
        
        M = [M ; vecM];
        
        uM = [uM;UvecM];
        
        k
    end
end
snpV;
%% start lda
close all
snpN = 1:numel(snpV);
initP = 6;
p = [];
h = [];
pv1 = [];
pv2 = [];
alpha = .05;
CL = {'r','g','b','r','g','b','r','g','b'};
for e = (initP+1):size(M,2)
    
    
    cnt = zeros(1,numel(snpN));
    for g = 1:numel(snpN)
        cnt(g) = sum(uM(:,e) == g);
    end
    
    
    [~,sidx] = sort(cnt,'descend');
    kp = sidx(1:2);
    f1 = uM(:,e) == kp(1);
    f2 = uM(:,e) == kp(2);
    
    
    subD1 = uM(f1,2:3);
    subD2 = uM(f2,2:3);
    group = [zeros(sum(f1),1);ones(sum(f2),1)];
    data = [subD1;subD2];
    
    
    lambda = myLDA(data,group);
    
    lambda = lambda/ norm(lambda);
    p1 = subD1*lambda;
    p2 = subD2*lambda;
    LAM(:,e-initP) = lambda;
    ANG(e-initP) = atan2(LAM(2,e-initP),LAM(1,e-initP));
    
    o1_1 = subD1(:,1);
    o1_2 = subD2(:,1);
    
     
    o2_1 = subD1(:,2);
    o2_2 = subD2(:,2);
    
    p(e-initP) = anova1([p1;p2],group,'off');
    pv1(e-initP) = anova1([o1_1;o1_2],group,'off');
    pv2(e-initP) = anova1([o2_1;o2_2],group,'off');
    
    %[h(e-initP),p(e-initP)] = ttest2(p1,p2);
    
    %[h1(e-initP),pv1(e-initP)] = ttest2(o1_1,o1_2);
    %[h2(e-initP),pv2(e-initP)] = ttest2(o2_1,o2_2);
     
    
    
    
    %{
    if mod(e,1000) == 0
        
        xx = 1:numel(p);
        for u = 1:9
            subplot(2,2,1);
            fidx = chr(1:numel(p)) == u;
            plot(xx(fidx),-log(p(fidx)),CL{u});hold on
            subplot(2,2,2);
            plot(xx(fidx),-log(pv1(fidx)),CL{u});hold on;
            subplot(2,2,3);
            plot(xx(fidx),-log(pv2(fidx)),CL{u});hold on;
            subplot(2,2,4);
            [fi,y] = ksdensity(ANG);
            plot(y,fi);
            axis([-pi pi 0 max(fi)])
            drawnow
            hold on
        end
        hold off
    end
    %}
    e
end

%%
[~,midx] = min(p);
e = midx;
    cnt = zeros(1,numel(snpN));
    for g = 1:numel(snpN)
        cnt(g) = sum(uM(:,e) == g);
    end
    
    
    [~,sidx] = sort(cnt,'descend');
    kp = sidx(1:2);
    f1 = uM(:,e) == kp(1);
    f2 = uM(:,e) == kp(2);
    
    
    subD1 = uM(f1,2:3);
    subD2 = uM(f2,2:3);
    group = [zeros(sum(f1),1);ones(sum(f2),1)];
    data = [subD1;subD2];
    
    
    lambda = myLDA(data,group);
    p1 = subD1*lambda;
    p2 = subD2*lambda;
    
    
    
%%
for e = 2:size(T,1)
    tmp = T{e,1};
    fidx = strfind(tmp,'_');
    key(e) = str2num(tmp((fidx(2)+1):end));
    e
end