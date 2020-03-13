FilePath = '/mnt/tetra/nate/stomataTopoData/Accessions_2016/';
FileList = {};
FileExt = {'nms'};
tic
FileList = gdig(FilePath,FileList,FileExt,1);
toc
%%
close all
I = imread(FileList{1000});
%I = rand(512,512);
% seed points
pt = randi(512,70,2);
pt(any(mod(pt,2)==1,2),:) = [];
imshow(I,[]);
hold on

pt = round(pt);
sim = zeros(size(I));



V = 0;
vTH = 20;
while V < vTH | min(D) <= 31
    dis = squareform(pdist(pt));
    diagIDX = sub2ind(size(dis),1:size(dis,1),1:size(dis,1));
    dis(diagIDX) = inf;
    [V,pidx] = min(dis(:));
    [q1,q2] = ind2sub(size(dis),pidx);
    
    D = abs(pt(q1,:) - pt(q2,:));
    
    
    if V < vTH | min(D) <= 31
        pt(q1,:) = [];
    end
    size(pt)
end



plot(pt(:,1),pt(:,2),'r*')
cp1 = pt(q2,:);
cp2 = pt(q1,:);
plot(cp1(1),cp1(2),'go');
plot(cp2(1),cp2(2),'mo');
box_size = abs(cp2 - cp1);
overlap_size = box_size/2;

overlap_box = [.5*(cp1 + cp2)-box_size/2 box_size];
box1 = [cp1-box_size 2*box_size];
box2 = [cp2-box_size 2*box_size];
delta_replace = [cp2-box_size - (.5*(cp1 + cp2)-box_size/2)];
replaceImage = imcrop(I,box2);
sim(box1(2):box1(2)+box1(4),box1(1):box1(1)+box1(3)) = imcrop(I,box1);

%%

randOverLapp = 10;
minOverLapp = 15;
close all
currentPoint = randi(512,1,2);
currentSize = randi(21,1,2) + 11;



    if currentSize(1) > currentPoint(1)
        currentSize(1) = currentPoint(1)+1;
    end
    
    
    if (512 - currentPoint(1)) < currentSize(1) 
        currentSize(1) = 512 - currentPoint(1);
    end
    
    
    if currentSize(2) > currentPoint(2)
        currentSize(2) = currentPoint(2)+1;
    end
    
    
    if (512 - currentPoint(2)) < currentSize(2) 
        currentSize(2) = 512 - currentPoint(2);
    end



sim = zeros(512,512);
simMask = zeros(512,512);
simMask(currentPoint(1)-currentSize(1):currentPoint(1)+currentSize(1),...
    currentPoint(2)-currentSize(2):currentPoint(2)+currentSize(2)) = 1;
imshow(simMask,[]);

seedImageNumber = randi(numel(FileList),1);
seedImageNumber = 1000;
I = imread(FileList{seedImageNumber});
imshow(I,[]);


% copy image from seed to sim
%box = [currentPoint-currentSize/2 2*currentSize+1];
%sim(box(2):box(2)+box(4),box(1):box(1)+box(3)) = imcrop(I,box);
sim(currentPoint(1)-currentSize(1):currentPoint(1)+currentSize(1),...
    currentPoint(2)-currentSize(2):currentPoint(2)+currentSize(2)) = ...
I(currentPoint(1)-currentSize(1):currentPoint(1)+currentSize(1),...
    currentPoint(2)-currentSize(2):currentPoint(2)+currentSize(2));
figure;
imshow(sim,[]);



%
for o = 1:35
FileList = FileList(randperm(numel(FileList)));
    % get next point
    DIST = bwdist(simMask);
    randD = randi(11,1,1) + 11;
    fidx = [];
    nextRegion = DIST > randD - 1 & DIST < randD + 1;
    [fidx(:,1),fidx(:,2)] = find(nextRegion==1);
    currentPointIDX = randi(size(fidx,1),1,1);
    currentPoint = fidx(currentPointIDX,:);

    currentSize = round(DIST(currentPoint(1),currentPoint(2)));
    currentSize = currentSize + minOverLapp + randi(randOverLapp,1,2);

    

    if currentSize(1) > currentPoint(1)
        currentSize(1) = currentPoint(1)+1;
    end
    
    
    if (512 - currentPoint(1)) < currentSize(1) 
        currentSize(1) = 512 - currentPoint(1);
    end
    
    
    if currentSize(2) > currentPoint(2)
        currentSize(2) = currentPoint(2)+1;
    end
    
    
    if (512 - currentPoint(2)) < currentSize(2) 
        currentSize(2) = 512 - currentPoint(2);
    end
    
    
    close all
    figure;
    imshow(cat(3,nextRegion,simMask,sim),[]);
    hold on
    plot(currentPoint(2),currentPoint(1),'y*');


    currentSearch = sim(currentPoint(1)-currentSize(1):currentPoint(1)+currentSize(1),...
        currentPoint(2)-currentSize(2):currentPoint(2)+currentSize(2));
    currentW = simMask(currentPoint(1)-currentSize(1):currentPoint(1)+currentSize(1),...
        currentPoint(2)-currentSize(2):currentPoint(2)+currentSize(2));
    figure;
    imshow(sim,[]);
    drawnow
    %%
    so = currentSearch;
    wo = currentW;



    sper = .25;
    fI = [];
    v = [];
    idx = [];
    parfor e = 1:500%1:100


            I = imread(FileList{e});
            oI = I;
            I = imresize(I,sper);
            s = imresize(so,sper);
            w = imresize(wo,sper);

            test = im2colF(I,size(s),[1 1]);
            delta = bsxfun(@minus,test,s(:));
            delta = bsxfun(@times,delta,w(:));
            delta = sum(delta.*delta,1);
            delta = col2im(delta,size(s),size(I),'sliding');
            [v(e),idx(e)] = min(delta(:));



            [q1,q2] = ind2sub(size(delta),idx(e));
            q1 = q1 * sper^-1;
            q2 = q2 * sper^-1;


            try
                fI(:,:,e) = oI(q2:q2+2*currentSize(1),...
                    q1:q1+2*currentSize(2));
            catch

            end
            e

    end



    [sv,sidx] = sort(v);
    toUse = NaN*zeros(size(sim));

    toUse(currentPoint(1)-currentSize(1):currentPoint(1)+currentSize(1),...
        currentPoint(2)-currentSize(2):currentPoint(2)+currentSize(2)) = fI(:,:,sidx(1));

    sim(sim==0) = NaN;
    sim = nanmean(cat(3,toUse,sim),3);
    imshow(sim,[]);
    sim(isnan(sim)) = 0;


    simMask(currentPoint(1)-currentSize(1):currentPoint(1)+currentSize(1),...
        currentPoint(2)-currentSize(2):currentPoint(2)+currentSize(2)) = 1;
end




%%


rectangle('Position',box1,'EdgeColor','r');
rectangle('Position',box2,'EdgeColor','b');
rectangle('Position',overlap_box,'EdgeColor','g');
so = imcrop(I,overlap_box);
so = imcrop(sim,overlap_box);
%s = zeros(size(s));
%s(4,:) = 1;
%s(:,4) = 1;
fI = [];
v = [];
sper = .25;
for e = 1:300%1:100
    
       
        I = imread(FileList{e});
        oI = I;
        I = imresize(I,sper);
        s = imresize(so,sper);
        
        
        test = im2colF(I,size(s),[1 1]);
        delta = bsxfun(@minus,test,s(:));
        delta = sum(delta.*delta,1);
        delta = col2im(delta,size(s),size(I),'sliding');
        [v(e),idx(e)] = min(delta(:));
        
        
        
        [q1,q2] = ind2sub(size(delta),idx(e));
        q1 = q1 * sper^-1;
        q2 = q2 * sper^-1;
        sbox = [q2 q1 box_size];


        rbox = [[q2 q1] + delta_replace 2*box_size];

        
        tmp = imcrop(oI,rbox);
        
        if all(size(tmp) == flip(2*box_size,2)+1)
        
        
            fI(:,:,e) = tmp;
            
            %{
            imshow(cat(2,fI(:,:,e),replaceImage),[]);

            drawnow
            %}
        end

        e
        
        
    
    
end



[sv,sidx] = sort(v);
toUse = NaN*zeros(size(sim));
toUse(box2(2):box2(2)+box2(4),box2(1):box2(1)+box2(3)) = fI(:,:,sidx(1));
sim(sim==0) = NaN;
sim = nanmean(cat(3,toUse,sim),3);



%%
close all
% seed point
pt = randi(512,1,2);
% seed box size
sz = randi(11,1,2) + 11;
npt = randi(512,1,2);
% copy data
sim(npt(1)-sz(1):npt(1)+sz(1),npt(2)-sz(2):npt(2)+sz(2)) = I(npt(1)-sz(1):npt(1)+sz(1),npt(2)-sz(2):npt(2)+sz(2));
% percent overlap
ov = .75;
% delta_point
d_nextpt(1) = randi(round(ov*sz(1)),1) + 5;
d_nextpt(2) = randi(round(ov*sz(2)),1) + 5;
% next point
pt = npt + [d_nextpt];
% next size
sz = randi(11,1,2) + 11;
box = [pt pt + 2*sz+1];
imshow(sim,[]);
rectangle('Position',box,'EdgeColor','r');
s = imcrop(sim,box);

