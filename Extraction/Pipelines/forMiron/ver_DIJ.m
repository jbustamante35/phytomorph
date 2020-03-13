%%
D = csvread('/home/nate/D1-1__lower=25_upper=200_smoothing=3__centroids.csv');
D = csvread('/home/nate/D1-2__lower=25_upper=200_smoothing=3__centroids.csv');
D = csvread('/home/nate/D1-3__lower=25_upper=200_smoothing=3__centroids.csv');
%%

close all
X = round(max(D(:,2)));
Y = round(max(D(:,3)));
M = zeros(X,Y);
for e = 1:size(D,1)
    M(round(D(e,3)),round(D(e,2))) = 1;
end
cM = imclose(M,strel('disk',21,0));
cM = imfill(cM,'holes');
cM = bwlarge(cM);
cM = imdilate(cM,strel('disk',31,0));
for e = 1:size(D,1)
    v(e) = cM(round(D(e,3)),round(D(e,2)));
end
fidx = find(v == 1);

figure;
imshow(cM,[])

dishT = 410;
dX = [];
uXY = mean(D(:,2:3));
for e = 1:size(D,1)
    dX(e,:) = D(e,2:3) - uXY;
end
delta = sum(dX.*dX,2).^.5;
figure;
plot3(D(:,2),D(:,3),D(:,1),'.')
%fidx = find(delta < dishT);
hold on
plot3(D(fidx,2),D(fidx,3),D(fidx,1),'r.')
view([0 90])
%%
D = D(fidx,:);
%%
close all
UQ = unique(D(:,1));
uA = [];
for e = 1:numel(UQ)
    fidx = find(D(:,1) == UQ(e));
    tmp = D(fidx,4:5);
    c1 = count(tmp(:,1));
    c2 = count(tmp(:,2));
    
    
    ct = (c1 == 1 & c2 == 1);
    if sum(ct) == 0
        ct = c1 == 1;
    end
    ct = c1 == 1;
    %ct = c1 == 2;
    
    
    
    uA(e) = mean(tmp(ct,1));
    uA(e) = median(tmp(ct,1));
    
    
    
    
    C(e) = sum(ct);
end

uC = mean(C);
mC = mode(C);

uAT = mean(uA);
uAT = median(uA);

% count with static area
for e = 1:numel(UQ)
    fidx = find(D(:,1) == UQ(e));
    tmp = D(fidx,4:5);
    
    tmpC = round(sum(tmp(:,1) .* uAT.^-1));
    tmpC = (sum(tmp(:,1) .* uAT.^-1));
    
    CT(e) = tmpC;
end

figure;
plot(C);
hold on;
plot(CT,'r');
title('Counts')

suA = imfilter(uA,fspecial('average',[1 21]),'replicate');
figure;
plot(uA)
hold on
plot(uAT*ones(size(uA)),'r')
plot(suA,'g')
title('Area')



% count with dynamic smoothed area
for e = 1:numel(UQ)
    fidx = find(D(:,1) == UQ(e));
    tmp = D(fidx,4:5);
    
    tmpC = round(sum(tmp(:,1) .* suA(e).^-1));
    tmpC = (sum(tmp(:,1) .* suA(e).^-1));
    
    CT2(e) = tmpC;
end



uCT = mean(CT);
uCT2 = mean(CT2);
[uC uCT uCT2 mean([uC uCT uCT2])]

%% dish check
figure;
plot3(D(:,2),D(:,3),D(:,1),'.')
view([0 90])
%%  tracer
%% scale the problem 
tmp = D(:,2:3);
[S C U E L ERR LAM] = PCA_FIT_FULL(tmp,2);
STD = diag(LAM).^-.5;
STD = STD';
C = bsxfun(@times,C,STD);
%%
D(:,2:3) = C;
%% backup D
Dbk = D;
%% restore D
D = Dbk;
%%
%% 
D = [D D(:,1)];
D(:,1) = .1*D(:,1);
%% look at some distances
UQ = unique(D(:,end));
for t = 1:(numel(UQ)-1)
    f1 = find(D(:,end) == UQ(t));
    f2 = find(D(:,end) == UQ(t+1));
    s1 = D(f1,1:3);
    s2 = D(f2,1:3);
    
    
    for e = 1:size(s1,1)
        delta = bsxfun(@minus,s2,s1(e,:));
        delta = sum(delta.*delta,2).^.5;
        mdSPACETIME{t}(e) = min(delta);
        
        delta = bsxfun(@minus,s2(:,2:3),s1(e,2:3));
        delta = sum(delta.*delta,2).^.5;
        mdSPACE{t}(e) = min(delta);
        
    end
    
    
    
    MdSPACETIME(t) = max(mdSPACETIME{t});
    MdSPACE(t) = max(mdSPACE{t});
    t
end
%% looking at min distance between objects helps this threshold for connections
R = 20;
A = Radjacency(D(:,1:3)',R);
%%
for e = 1:numel(UQ)
    fidx = find(D(:,end) == UQ(e));
    for i = 1:numel(fidx)
        for j = 1:numel(fidx)
             A(fidx(i),fidx(j)) = 0;
        end
    end
    e
end
%%
close all
fidxStr = find(D(:,end) == 0);
fidxStp = find(D(:,end) == max(UQ));

for str = 1:numel(fidxStr)
    for stp = 1:numel(fidxStp)
        
     
        
        
        [path{str,stp} , pathcost(str,stp)] = dijkstra(A , fidxStr(str) , fidxStp(stp));
       
        
        plot(D(:,2),D(:,3),'b.','MarkerSize',1)
        hold on
        plot(D(fidxStr(str),2),D(fidxStr(str),3),'go','MarkerSize',5)
        plot(D(fidxStp(stp),2),D(fidxStp(stp),3),'ro','MarkerSize',5)
        plot(D(path{str,stp} ,2),D(path{str,stp} ,3),'m','LineWidth',3)
        hold off
        drawnow
        pause(.4)
        
        
        
    end
end







