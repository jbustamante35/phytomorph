%% find the HSV video files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FilePath = '/mnt/spaldingdata/nate/octerineDataStore/planarian/';
FileList = {};
FileExt = {'hsv'};
tic
FileList = fdig(FilePath,FileList,FileExt,1);
toc
%% read counts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cData = readtext('/mnt/spaldingdata/nate/octerineDataStore/planarian/counts.csv');
header = cData(:,1);
timeHeader = cell2mat(cData(1,2:(end-1)));
timeCounts = cell2mat(cData(2:end,2:(end-1)));
totalCounts = cell2mat(cData(2:end,end));
%% find the main dish
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for e = 1:numel(FileList)
    videoFile = FileList{e};
    slice = readSlice(videoFile,1,0);
    imshow(slice,[]);
    drawnow
    waitforbuttonpress
end
%% get the petri masks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
skipN = 100;
sigSmoothValue = 11;
per = .4;
edgeCloseValue = 31;
edgeDilateValue = 3;
disp = false;
parfor e = 1:numel(FileList)
    videoFile = FileList{e};
    [maskStack(:,:,e),experimentFrame(e),maskStruct(e),background(:,:,e)] = findDishMask(videoFile,skipN,sigSmoothValue,per,edgeCloseValue,edgeDilateValue,disp);
     if disp
    %    out = flattenMaskOverlay(maskStruct(e).testImage,maskStruct(e).maskImage,.4,'r');
    %    imshow(out,[]);
    %    drawnow
    end
    fprintf(['done with:' num2str(e) ':' num2str(numel(FileList)) '\n']);
end
%{
%% look at temporal signal in dish over time
% not much trend over time
% rather I fixed the background code to not include the true zero signal
figure;
close all
e = 1;
videoFile = FileList{e};
[tmpSlice,szC] = readSlice(videoFile,1,0);

str = experimentFrame(e);
stp = szC(3);

fidx = find(maskStack(:,:,e)==1);
cnt = 1;
skipN = 5;
sig = [];
for t = str:skipN:stp
    [tmpSlice,szC] = readSlice(videoFile,t,0);
    sig(cnt) = mean(tmpSlice(fidx));
    cnt = cnt + 1;
    plot(sig)
    drawnow
end
%}
%% loop over frames and view petri plate mask
% date: March, 30 2020
% reviewing petri plate masks with josh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for e = 1:size(maskStack,3)
    %videoFile = FileList{e};
    %tmpSlice = readSlice(videoFile,100,0);
    %tmpSlice = squeeze(tmpSlice);
    tmpSlice = background(:,:,e);
    tmpMask = maskStack(:,:,e) > .4;
    dB = bwboundaries(tmpMask);
    out = flattenMaskOverlay(tmpSlice,tmpMask,.1,'r');
    imshow(out,[]);hold on
    plot(dB{1}(:,2),dB{1}(:,1),'b')
    drawnow
    hold off
    waitforbuttonpress
end
%{
%% test the optical flow from matlab
% testing worked
% date: March 31,2020
e = 1;
opticFlow1 = opticalFlowFarneback('NeighborhoodSize',15);
videoFile = FileList{e};
t = 0;
c = [];
r = [];
for t = 0:1:50
   
    I = readSlice(videoFile,experimentFrame(e)+t,0);

    if t == 0
        [c(t+1),r(t+1),V] = impixel(I);
    end

    spacing = 10;
    if t == 0
        [g1,g2] = ndgrid(1:size(I,1),1:size(I,2));
        grid = (mod(g1,spacing) == 0) & (mod(g2,spacing) == 0);
        fidx = find(grid & maskStack(:,:,e));
        x = g2(fidx);
        y = g1(fidx);
    end


    flowField = estimateFlow(opticFlow1,I);


    u = flowField.Vx(fidx);
    v = flowField.Vy(fidx);

    c(t+2) = c(t+1) + ba_interp2(flowField.Vx,c(t+1),r(t+1));
    r(t+2) = r(t+1) + ba_interp2(flowField.Vy,c(t+1),r(t+1));

    imshow(I,[]);
    hold on
    quiver(x,y,u,v,5)
    plot(c(end),r(end),'r.')
    drawnow
    hold off
end
%}
%% gather the stats on velocity intensity
clear loaderFunc
loaderFunc{1} = @(x,ni)randomVIslice(x,5,15);
loaderFunc{2} = @(x,ni,I)find(maskStack(:,:,ni)==1);

I = loaderFunc{1}(FileList{1});
dX = {linspace(0,1,255),linspace(0,15,70)};
[density] = channelSampler(FileList,50,10000,dX,loaderFunc);
%% apply to image
close all
I = loaderFunc{1}(FileList{1});
pI = density.prob(I);
imshow(pI,[])
%% 
imshow(density.density,[])

%% test the optical flow from matlab - main here

% date: March 31,2020
e = 1;


velocityParticleFilter = 75;
velocityCloseValue = 5;


objectParticleFilter = 50;
objectParticleFilterBEST = 90; % best is moving by not object-detected added to object-detected
objectParticleFilterBEST = 50;

N = {};
V = {};
modN = 1;
disp = true;
collectFeatures = true;
ON = [];
VN = [];
ON_BEST = [];
close all
I = background(:,:,1);


gaugeF = zeros(size(I));
gaugeF(378,118) = 1;
gaugeF = imdilate(gaugeF,strel('disk',51,0));
gidx = find(gaugeF);
gB = bwboundaries(gaugeF);

samW = 30;
samN = 15;
[sam1,sam2] = ndgrid(linspace(-samW,samW,2*samN+1),linspace(-samW,samW,2*samN+1));
samD = [sam1(:),sam2(:),ones(size(sam1(:)))];
domain.x = samD;
domain.sz = size(sam1);

objectStore = {};

clear staicStore dynamicStore;
cntO = 1;

for e = 1:10%1:numel(FileList)


    tmC = timeCounts(e,:);
    %tmC = interp1(1:numel(tmC),tmC,linspace(1,numel(tmC),1800),'cubic');


    opticFlow1 = opticalFlowFarneback('NeighborhoodSize',15);
    %opticFlow1 = opticalFlowLK();
    %opticFlow1 = opticalFlowLKDoG();

    videoFile = FileList{e};


    [~,M] = readSlice(videoFile,1,0);
    M(3) = 1800;
    totalFrames = M(3) - experimentFrame(e);
    
    
    totalFrames = 50;
    
    
    velFrame = [];
    particleNumber = [];
    particleSpeed = {};
    velStack = [];
    iStack = [];
    gv = [];

    ON = [];
    ON_BEST = [];
    VN = [];
    tmpObjectStore = {};
    
    
    
    for t = 1:totalFrames
        try
            frameInd = experimentFrame(e) + t;
            tic;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % read the image(s)
            I = readSlice(videoFile,experimentFrame(e) + t,0);
            % get the petri mask from first pass methed(s)
            petriMask = logical(maskStack(:,:,e));
            % get the bounding box for the petri plate
            pR = regionprops(petriMask);
            % crop mask and image and background
            I = imcrop(I,pR.BoundingBox);
            % crop the petri mask
            petriMask = logical(imcrop(petriMask,pR.BoundingBox));
            % crop the backgrond
            bg = imcrop(background(:,:,e),pR.BoundingBox);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % object detection
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % remove the background
            objectImage = I - bg;
            % image - background
            I_lessBK = objectImage;
            % calculate the threshold
            thresh = graythresh(objectImage(petriMask));
            % threshold the image
            objects = objectImage(petriMask) > thresh;
            % create object mask
            objectImage = zeros(size(objectImage));
            % pouring threshold results
            objectImage(petriMask) = objects;

            %%%%%%%%%%%%%%%
            % local does not seem to help
            % rather it hurts
            %%%%%%%%%%%%%%%
            localT = false;
            if localT
                objectImage = bindVec(I - bg);
                thresh = adaptthresh(objectImage);
                objectImage = imbinarize(objectImage,thresh);
                objectImage = objectImage .* petriMask;
                %imshow(objectImage,[]);
                %drawnow;waitforbuttonpress;
            end
            %%%%%%%%%%%%%%%
            
           
            baseLineObjectImage = logical(objectImage);
            % small area filter the object
            objectImage = bwareaopen(objectImage,objectParticleFilter);
            % calc distance map to objects
            objectDistance = double(bwdist(objectImage));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % flow data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % velocity field
            flowField = estimateFlow(opticFlow1,I);
            % calc flow field
            u = flowField.Vx;
            v = flowField.Vy;
            % calc the speed
            s = (u.^2 + v.^2).^.5;
            % calc the diergence
            div = divergence(u,v);
            % calc the integral of div over gauge area
            gv(t) = mean(div(gidx));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % sample the image (gray scale) and speed
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            in = I(petriMask);
            sv = s(petriMask);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % make velocity mask
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % normalize the velocity
            svn = bindVec(sv);
            vm = svn > graythresh(svn);
            vMask = zeros(size(petriMask));
            vMask(petriMask) = vm;
            vMask = logical(vMask);
            % close mask for veloctiy
            vMask = imclose(vMask,strel('disk',velocityCloseValue,0));
            % size filter on velocities
            vMask = bwareaopen(vMask,velocityParticleFilter);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get the data from the binary masks
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            velR = regionprops(vMask,'PixelIdxList','Area','Centroid','BoundingBox','Orientation');
            objR = regionprops(logical(objectImage),'PixelIdxList','Area','Centroid','BoundingBox','Orientation');
            
            %comMask = logical(objectImage) | vMask;
            %comR = regionprops(comMask,'PixelIdxList','Area','Centroid','BoundingBox','Orientation');
            
            % counts from object and velocity blobs
            ON(t) = numel(objR);
            VN(t) = numel(velR);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % START my track
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
break
%%
close all
            domainW = 7;
            domainN = 2*domainW + 1;
            X = linspace(-domainW,domainW,domainN);
            Y = linspace(-domainW,domainW,domainN);
            x = basisT.affineSpace2D(X,Y);
            
            X = linspace(0,domainW,domainW+1);
            Y = linspace(-pi,pi,50);
            %x = basisT.affineSpace2D_disk(X,Y);

            opticFlow1 = opticalFlowFarneback('NeighborhoodSize',15);
            prev = [];
            
           	e = 1;t = 1;clear T TT;
            for t = 1:4%30


                zp = basisT([0;0;1]);

                
                I = basisF(readSlice(videoFile,experimentFrame(e) + t,0));
                J = basisF(readSlice(videoFile,experimentFrame(e) + t+1,0));
                
           
                
                if t == 1
                    K = basisF(readSlice(videoFile,experimentFrame(e),0));
                    flowField = estimateFlow(opticFlow1,I.E);
                    %{
                    for h = 1:30
                        I = basisF(readSlice(videoFile,experimentFrame(e) + h,0));
                        imshow(I)
                    end
                    %}
                    %[c,r,V] = impixel(I.E);
                    T(t) = basisT.affine2D([c;r]);
                    TT(t) = basisT.affine2D([c;r]);
                end
                flowField = estimateFlow(opticFlow1,J.E);
                
                zoom_pt = T(t)*zp;
                
                
                boxW = 100;
                box = zoom_pt.point2box([boxW boxW]);
                i1 = I*box;
                i2 = J*box;
                
                %{
                testPattern1 = zeros(size(i1.E));
                testPattern1(51,51) = 1;
                testPattern2 = zeros(size(i1.E));
                testPattern2(55,51) = 1;
                i1 = basisF(double(bwdist(testPattern1)));
                i2 = basisF(double(bwdist(testPattern2)));
                %}
                
                T1 = i1.generateCenterAffine();
                T2 = i2.generateCenterAffine();

                %{
                sub1 = i1*(T1*x);
                sub2 = i2*(T2*x);
                imshow(i1);
                hold on;
                hsz = hsize(i1);
                plot(hsz(2),hsz(1),'r*');
                figure;imshow(sub1);
                %}
                
                
                % generate function to generate transformation
                dX_range = [[20;-20],[-20;20]];
                dR_range = [pi;-pi]/2;
                dS_range = [[1.7;.5],[1.7;.5]];
                [K,initX] = basisT.makeTF(dX_range,dR_range,dS_range,[0 4 0 1 1]);
                [K,initX,mf] = basisT.makeTF(dX_range,dR_range,dS_range);
                
                
                %P = K(initX);
                %P.E
                
                K = @(v)basisT.transform2d2(v);
                initX = [0 0 pi/4 1 1 pi/8];
                
                %if ~isempty(prev);initX = prev;end
                
                sub1 = prod([i1 T1 x]);
                sub2 = @(v)prod([i2 T2 K(v) x]);
                %{
                K = @(v)basisT(reshape(v,[3 3]));
                sub2 = @(v)prod([i2 T2 K(v) x]);
                initX = eye(3);
                initX = initX(:)';
                %}
                
                
                
                contrast = @(v)norm(sub1 - sub2(v));
                
                %{
                m1 = [K(xSol) T2 sp];
                v1 = prod(m1(1:end),2);
                sp = basisT([10;10;1]);
                P = K(xSol)*(T2*sp);
                i2*P
                v1.E;
                P.E
                
                qp = basisT([3;3;1]);
                w = K(xSol);
                where = basisT(inv(w.E))*qp;
                where.E
                %}
                solverType = 'patternsearch';
                solverType = 'fminsearch';
                solverType = 's';
                jumpMAX = 15;
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
                
                
                %mf(xSol)
                xSol
                SOL = K(xSol);
                SOL.E
                SOL.E(end,:) = [0 0 1];
                prev = xSol;
                %SOL.E = inv(SOL.E);
                %P.E = inv(P.E);
                
                orn  = [cos(xSol(3)) sin(xSol(3))];
                
                T(t+1) = T(t)*SOL;
                
                
                dx1 = ba_interp2(flowField.Vx,TT(t).E(1,3),TT(t).E(2,3));
                dx2 = ba_interp2(flowField.Vy,TT(t).E(1,3),TT(t).E(2,3));
                PP = eye(3);
                PP(1,3) = dx1;
                PP(2,3) = dx2;
                
                
                %[xSol,v3] = patternsearch(contrast, [dx1 dx2 0 1 1],A,b,Aeq,beq,lb,ub,nonlcon,options);
                        
                
                
                
                SOL2 = basisT(PP);
                TT(t+1) = TT(t)*SOL2;
                
                
                sp = basisT([10;10;1]);
                q = K(xSol)*(T2*sp);
                q.E;
                
                toP = T(t)*zp;
                toPlot = T(t+1)*zp;
                rec = toPlot.point2box([size(sub1.E,1) size(sub1.E,1)]);
                recO = toP.point2box([size(sub1.E,1) size(sub1.E,1)]);
                toPlot2 = TT(t+1)*zp;
                RGB = cat(3,I.E,J.E,zeros(size(I.E)));
                
                
                imshow(J.E,[]);hold on
                plot(toPlot.E(1),toPlot.E(2),'r.');
                plot(toPlot2.E(1),toPlot2.E(2),'g.');
                plot(toP.E(1),toP.E(2),'b.');
                quiver(toP.E(1),toP.E(2),orn(1),orn(2),80,'Color','b')
                drawnow
                %
                
                test = prod([i2 T2 x]);
                A = sub1.E;
                B = sub2(xSol);
                B = B.E;
                C = test.E;
                imshow(interp2([A,B,C],2),[])
                %imshow(([A,B,C]),[])
                %{
                rectangle('Position',rec.E,'EdgeColor','r');
                rectangle('Position',recO.E,'EdgeColor','b');
                %}
                drawnow
               
                %waitforbuttonpress
                %{
                test = prod([i2 T2 x]);
                A = sub1.E;
                B = sub2(xSol);
                B = B.E;
                C = test.E;
                figure;imshow(A);
                figure;imshow(B);
                figure;imshow(C);
                waitforbuttonpress;
                close all
                %}
            end
%%
            t = prod([i2 T2 x]);
            j1 = sub1;
            j2 = sub2(xo);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % END my track
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %}   
 

            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % inline display code - not needed anymore
            mag = 30;
            imshow(I,[]);hold on
            for o = 1:numel(objR)
                TAN = mag*[sin(objR(o).Orientation*pi/180);cos(objR(o).Orientation*pi/180)];
                NOR = [TAN(2);-TAN(1)];
                DX = [objR(o).Centroid]';
                F = [TAN NOR DX];
                plot(DX(1),DX(2),'g.')
                quiver(F(1,3),F(2,3),F(1,1),F(2,1),'r')
                quiver(F(1,3),F(2,3),F(1,2),F(2,2),'b')
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %}
          
            if collectFeatures
                [velR] = objectFeatures(I_lessBK,cat(3,u,v),objectImage,vMask,velR,domain);
                [objR] = objectFeatures(I_lessBK,cat(3,u,v),objectImage,vMask,objR,domain);
                %[comR] = objectFeatures(I_lessBK,cat(3,u,v),objectImage,vMask,comR,domain);
                %tmpObjectStore{t}.velocityBlobs = velR;
                %tmpObjectStore{t}.objectBlobs = objR;
            end
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % work up on velocity object
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for vo = 1:numel(velR)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % sample the centroid velocity
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                velR(vo).cu = ba_interp2(u,velR(vo).Centroid(1),velR(vo).Centroid(2));
                velR(vo).cv = ba_interp2(v,velR(vo).Centroid(1),velR(vo).Centroid(2));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % make a velocity reference frame
                % and measure
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                TAN = [velR(vo).cu velR(vo).cv];
                TAN = TAN / norm(TAN);
                NOR = [TAN(2) -TAN(1)];
                % get the vel field in the blob
                tmpVelField = [u(velR(vo).PixelIdxList),v(velR(vo).PixelIdxList)];
                % get velocities relative to the frame
                alongDir = tmpVelField*TAN';
                alongNor = tmpVelField*NOR';
                % measure the expected velocities - (x,y)
                velR(vo).du = mean(u(velR(vo).PixelIdxList));
                velR(vo).dv = mean(v(velR(vo).PixelIdxList));
                % measure the expected velocities - (u,v)
                velR(vo).alongT = mean(alongDir);
                velR(vo).alongN = mean(alongNor);
                % measure the std velocities - (u,v)
                velR(vo).alongT_var = std(alongDir);
                velR(vo).alongN_var = std(alongNor);
                % measure the centroid
                velR(vo).x = velR(vo).Centroid(1);
                velR(vo).y = velR(vo).Centroid(2);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
            
            if collectFeatures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % for plotting and display
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                cenV = [[velR.centroid_vx]',[velR.centroid_vy]'];
                avgV = [[velR.expected_vx]',[velR.expected_vy]'];
                cenSpeed = sum(cenV.*cenV,2).^.5;
                avgSpeed = sum(avgV.*avgV,2).^.5;
                corSpeed = sum(cenV.*avgV,2).*(avgSpeed.^-1).*(cenSpeed.^-1);
                alongT = [velR.expected_vt]';
                alongN = [velR.expected_vn]';
                pos = [[velR.x]',[velR.y]'];
                varV = [[velR.var_vt]',[velR.var_vn]'];
                cir = [cos(linspace(-pi,pi))',sin(linspace(-pi,pi))'];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            plot(avgSpeed);hold on
            plot(cenSpeed);
            plot(corSpeed);
            plot(alongT);
            plot(alongN);
            hold off
            legend({'average','centroid','correlation','alongT','alongN'});
            pause(.5);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % find moving objects that are not found via background sub
            % optical flow finds objects but background subtraction does
            % not
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            toObject = [];
            socialDistance = 50;
            for i = 1:numel(velR)
                toObject(i) = ba_interp2(objectDistance,velR(i).Centroid(1),velR(i).Centroid(2));
            end
            % if object is "too" far away from intensity object
            oidx = find(toObject > socialDistance);
            z1 = zeros(size(objectImage));
            for c = 1:numel(oidx)
                z1(velR(c).PixelIdxList) = 1;
            end
            % add velocity objects to the intensity
            bestObjectImage = objectImage | z1;
            bestObjectImage = bwareaopen(bestObjectImage,objectParticleFilterBEST);
            BESTobjR = regionprops(logical(bestObjectImage),'PixelIdxList','Area','Centroid');
            ON_BEST(t) = numel(BESTobjR);
            
            comR = regionprops(logical(bestObjectImage),'PixelIdxList','Area','Centroid','BoundingBox','Orientation');
            
            [comR] = objectFeatures(I_lessBK,cat(3,u,v),objectImage,vMask,comR,domain);
               
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %{
            %%%%%%%%%%%%%%%%%%%%%
            % tried counting - didn't work
            %%%%%%%%%%%%%%%%%%%%%
            cidx = count([objR.Area]);
            fidx1 = find(cidx==1);
            fidx2 = find(cidx==2);
            z1 = zeros(size(objectImage));
            z2 = zeros(size(objectImage));
            for c = 1:numel(fidx1)
                z1(objR(c).PixelIdxList) = 1;
            end
            for c = 1:numel(fidx2)
                z2(objR(c).PixelIdxList) = 1;
            end
            %}
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % store velocities
            velStack(:,t) = sv(:);
            % store the intensities
            iStack(:,t) = in(:);
            %
            particleSpeed{t} = [];
            particleNumber(t) = numel(velR);
            for r = 1:numel(velR)
                particleSpeed{t}(r) = mean(s(velR(r).PixelIdxList));
            end
            velFrame(t) = mean(particleSpeed{t});
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if disp

                
                
                if mod(t,modN) == 0
                    out = I;
                    out = flattenMaskOverlay(out,baseLineObjectImage,.5,'r');
                    out = flattenMaskOverlay(out,logical(objectImage),.5,'g');
                    out = flattenMaskOverlay(out,logical(bestObjectImage),.8,'b');

                    out = flattenMaskOverlay(out,logical(vMask),.5,'c');

                    %out = flattenMaskOverlay(out,logical(z1),.5,'c');
                    %out = flattenMaskOverlay(out,logical(z2),.5,'y');
                    out = flattenMaskOverlay(out,logical(petriMask),.1,'y');
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % work on up display for finding single/isolated schistosome
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % collect data for up-front processing for obj NOT vel
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % for plotting and display
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                toOp = comR;
                cenV = [[toOp.centroid_vx]',[toOp.centroid_vy]'];
                avgV = [[toOp.expected_vx]',[toOp.expected_vy]'];
                cenSpeed = sum(cenV.*cenV,2).^.5;
                avgSpeed = sum(avgV.*avgV,2).^.5;
                corSpeed = sum(cenV.*avgV,2).*(avgSpeed.^-1).*(cenSpeed.^-1);
                alongT = [toOp.expected_vt]';
                alongN = [toOp.expected_vn]';
                pos = [[toOp.x]',[toOp.y]'];
                varV = [[toOp.var_vt]',[toOp.var_vn]'];
                cir = [cos(linspace(-pi,pi))',sin(linspace(-pi,pi))'];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                deltaVec = [];deltaDis = [];deltaMat = [];
                for source = 1:size(pos,1)
                    parfor target = 1:size(pos,1)
                        deltaVec(source,target,:) = pos(target,:) - pos(source,:);
                        deltaDis(source,target) = norm(squeeze(deltaVec(source,target,:)));
                        deltaMat(source,target,:,:) = [pos(target,:);pos(source,:)];
                    end
                end
                [sortedDistance,sidx] = sort(deltaDis,2);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % work on up display for finding single/isolated schistosome
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %subplot(1,2,1);
                
                
                imshow(out,[]);
                drawnow
                title(num2str(t));
                hold on
                plot(pos(:,1),pos(:,2),'r.')
                %{
                for e = 1:size(pos,1)
                    text(pos(e,1),pos(e,2),num2str(e),'BackgroundColor','w')
                end
                %}
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % threshold on velocity and distance 
                % currently this will sample movers
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                distanceT = 40;
                velocityT = 4;
                bulkV = 20;
                for source = 1:size(pos,1)
                    if cenSpeed(source) < velocityT
                        %for target = 1:size(pos,1)
                            if all(sortedDistance(source,2:end) > distanceT)
                                
                                bulkBox = toOp(source).BoundingBox;
                                bulkBox(1:2) = bulkBox(1:2) - bulkV;
                                bulkBox(3:4) = bulkBox(3:4) + 2*bulkV;
                                
                                % crop 
                                l1 = imcrop(I,bulkBox);
                                l2 = imcrop(I_lessBK,bulkBox);
                                l3 = imcrop(objectImage,bulkBox);
                                
                                
                                staticStore(cntO).I = cat(3,l1,l2,l3);
                                staticStore(cntO).P = [pos(source,:) frameInd];
                                cntO = cntO + 1;
                                
                                
                                rectangle('Position',bulkBox,'EdgeColor','r');
                                % plot the two shortest = three shortest - self
                                for n = 2:3
                                    plot(squeeze(deltaMat(source,sidx(source,n),:,1)),squeeze(deltaMat(source,sidx(source,n),:,2)),'c')

                                    %plot(squeeze(deltaMat(source,target,:,1)),squeeze(deltaMat(source,target,:,2)),'c')

                                end
                            end
                        %end
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % threshold on velocity and distance 
                % currently this will sample movers
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                distanceT = 40;
                velocityT = 4;
                bulkV = 20;
                for source = 1:size(pos,1)
                    if cenSpeed(source) > velocityT
                        %for target = 1:size(pos,1)
                            if all(sortedDistance(source,2:end) > distanceT)
                                bulkBox = toOp(source).BoundingBox;
                                bulkBox(1:2) = bulkBox(1:2) - bulkV;
                                bulkBox(3:4) = bulkBox(3:4) + 2*bulkV;
                                
                                % crop 
                                l1 = imcrop(I,bulkBox);
                                l2 = imcrop(I_lessBK,bulkBox);
                                l3 = imcrop(objectImage,bulkBox);
                                
                                % 
                                dynamicStore(cntO).I = cat(3,l1,l2,l3);
                                dynamicStore(cntO).P = [pos(source,:) frameInd];
                                cntO = cntO + 1;
                                
                                % 
                                rectangle('Position',bulkBox,'EdgeColor','g');
                                % plot the two shortest = three shortest - self
                                for n = 2:3
                                    plot(squeeze(deltaMat(source,sidx(source,n),:,1)),squeeze(deltaMat(source,sidx(source,n),:,2)),'c')

                                    %plot(squeeze(deltaMat(source,target,:,1)),squeeze(deltaMat(source,target,:,2)),'c')

                                end
                            end
                        %end
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                
                
               
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % line for quiver velocity
                %quiver(pos(:,1),pos(:,2),cenV(:,1),cenV(:,2),'r')
                mag = 5;
                for q = 1:size(pos,1)
                    T = cenV(q,:);
                    T = T / norm(T);
                    NOR = [T(2) -T(1)];
                    F = mag*[T' NOR'];
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%
                    % make yellow circle for variance
                    ctmp(:,1) =  cir(:,1)*varV(q,1);
                    ctmp(:,2) =  cir(:,2)*varV(q,2);
                    ctmp = (F*ctmp')';
                    %plot(ctmp(:,1)+pos(q,1),ctmp(:,2)+pos(q,2),'y')
                    %%%%%%%%%%%%%%%%%%%%%%%%
                end
                %%% gauge area
                %plot(gB{1}(:,2),gB{1}(:,1),'g')
                %%% gauge area
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                hold off
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% line graph start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               
                %{
                subplot(1,2,2);
                % gauge plots
                igv = cumsum(gv);
                plot(igv);
                hold on
                plot(t,igv(t),'r*')
                % gauge plots
                %}

                %{
                subplot(1,2,2)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                tvec = (1:t) - 1;
                %plot(tvec,ON,'k');hold on
                %plot(tvec,VN,'r');
                %plot(tvec,ON_BEST,'c');



                av = .5*mean(ON_BEST)+.5*mean(ON);
    



              

                plot(tvec,totalCounts(e)*ones(size(ON)),'k','LineWidth',1);

                hold on
                plot(tvec,mean(ON_BEST)*ones(size(ON)),'c--')
                plot(tvec,mean(ON)*ones(size(ON)),'k--')
                plot(tvec,(av)*ones(size(ON)),'g--');
               



                smoothAmount = 11;
                notMoving1 = .5*(ON_BEST+ON) - VN;
                %plot(tvec,notMoving1,'b-','LineWidth',.5)
                smoothNotMoving1 = imfilter(notMoving1,ones(1,smoothAmount)/smoothAmount,'replicate');
                plot(tvec,smoothNotMoving1,'b','LineWidth',1);


                notMoving2 = ON_BEST - VN;
                %plot(tvec,notMoving2,'c-','LineWidth',.5)
                smoothNotMoving2 = imfilter(notMoving2,ones(1,smoothAmount)/smoothAmount,'replicate');
                plot(tvec,smoothNotMoving2,'c','LineWidth',1);

                %{
                notMoving3 = mean(ON_BEST)*ones(size(ON_BEST)) - VN;
                %plot(tvec,notMoving2,'c-','LineWidth',.5)
                smoothNotMoving3 = imfilter(notMoving3,ones(1,smoothAmount)/smoothAmount,'replicate');
                plot(tvec,smoothNotMoving2,'c','LineWidth',1);
                %}


                smoothAmount = 5;
                smoothV = imfilter(VN,ones(1,smoothAmount)/smoothAmount,'replicate');
                notMoving3 = av*ones(size(ON_BEST)) - smoothV;
                %plot(tvec,notMoving2,'c-','LineWidth',.5)
                smoothNotMoving3 = imfilter(notMoving3,ones(1,smoothAmount)/smoothAmount,'replicate');
                plot(tvec,smoothNotMoving3,'m','LineWidth',1);




                sam = 1:50:numel(tvec);
                %sam = timeHeader*10;
                plot(tvec(sam),smoothNotMoving2(sam),'r')
                plot(tvec(sam),smoothNotMoving2(sam),'ro')

                plot(timeHeader*10,tmC,'k','LineWidth',1);
                plot(timeHeader*10,tmC,'k*');
        
                hold off
                axis([0 t -10 50]);
                %}
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% line graph end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                hold off
                drawnow

                %{
                for o = 1:numel(velR)
                    text(objR(o).Centroid(1),objR(o).Centroid(2),num2str(objR(o).Area),'BackgroundColor','w');
                end
                %}
                %{
                subplot(1,2,1);
                imshow(out,[]);
                drawnow

                subplot(1,2,2)
                yyaxis('left')
                plot(velFrame(min(numel(velFrame),2):end))
                yyaxis('right')
                plot(particleNumber(min(numel(velFrame),2):end))
                title(num2str(t))
                drawnow
                %}
            end
            fprintf(['done frame:' num2str(t) ':' num2str(totalFrames) ':' num2str(e) '\n'])
        catch ME
            getReport(ME)
        end
    end



    objectStore{e} = tmpObjectStore;

    N{e} = particleNumber;
    V{e} = velFrame;

    NUM_BEST{e} = ON_BEST;
    NUM_IN{e} = ON;
    NUM_VEL{e} = VN;

    

end
%%
szS = [];
for e = 1:numel(staticStore)
    if ~isempty(staticStore(e).I)
        szS(e,:) = size(staticStore(e).I);
    else
        szS(e,:) = [0 0 0];
    end
end
rmidx = all(szS == 0,2);
staticStore(rmidx) = [];
%%
[X,Y,Z,A] = synthesizeImage(background,maskStack,staticStore,30,200);
%%
poolobj = gcp('nocreate');
delete(poolobj);
%%
close all force
per = .25;
smSZ = 21;
flt = zeros(smSZ);
flt((numel(flt)+1)/2) = 1;
flt = bwdist(flt);
flt = flt <= (smSZ-1)/2;
flt((numel(flt)-1)/2) = 1;

            subX = [];
            trainY = Y;

            newX = [];
            newY = [];
            
            for e = 1:size(X,3)
                tmp = imresize(X(:,:,e),per);
                tQ = im2colF(tmp,[smSZ smSZ],[1 1]);
                tQ = bsxfun(@times,flt(:),tQ);
                tQ = sort(tQ,1);
               
                
                tmp = imresize(double(edge(X(:,:,e))),per);
                tQ2 = im2colF(tmp,[smSZ smSZ],[1 1]);
                tQ2 = bsxfun(@times,flt(:),tQ2);
                tQ2 = sort(tQ2,1);
                
                
                tmp = imresize(A(:,:,e),per);
                tP = im2colF(tmp,[21 21],[1 1]);
                tP = abs(sum(tP,1));
                
                
                tP = discretize(tP,0:1:40);
                
                UQ = unique(tP);
                TOT = [];
                for v = 1:numel(UQ)
                    TOT(v) = sum(tP==UQ(v));
                    
                    fidx = find(tP == UQ(v));
                    fidx = fidx(randperm(numel(fidx)));
                    toS = min(numel(fidx),20);
                    tX = [tQ(:,fidx(1:toS));tQ(:,fidx(1:toS))];
                    newX = [newX tX];
                    newY = [newY UQ(v)*ones(1,toS)];
                end
                size(newX)
                
            end
            %%
            
            [S C U E L ERR LAM] = PCA_FIT_FULL_T(newX,5);
            cNet = feedforwardnet([2 18]);
            cNet = train(cNet,C,newY);
            subX2 = [];
            for e = 1:size(X,3)
                tmp = imresize(X(:,:,e),per);
                tmp2 = imresize(double(edge(X(:,:,e))),per);
                
                tQ = im2colF(tmp,[smSZ smSZ],[1 1]);
                tQ = bsxfun(@times,flt(:),tQ);
                tQ = sort(tQ,1);
                
                tQ2 = im2colF(tmp2,[smSZ smSZ],[1 1]);
                tQ2 = bsxfun(@times,flt(:),tQ2);
                tQ2 = sort(tQ2,1);
                
                
                tQ = PCA_REPROJ_T([tQ;tQ2],E,U);
                tQ = sim(cNet,tQ);
                tQ = col2im(tQ,[21 21],size(tmp));
                %imshow(tQ,[]);
                %drawnow
                subX2(:,:,e) = tQ;
                e
            end

            
           %% 
close all force
per = .25;
            subX = [];
            trainY = Y;

            newX = [];
            newY = [];
            
            for e = 1:size(X,3)
                subX(:,:,e) = imresize(X(:,:,e),per);
            end

            
            szX = size(subX);
            subX = reshape(subX,[szX(1) szX(2) 1 szX(3)]);

            
            szX = size(subX2);
            subX = reshape(subX2,[szX(1) szX(2) 1 szX(3)]);

            
            [sY,sU,sS] = zscore(Y);
            cv = cvpartition(size(subX,4),'HoldOut',.2);
            trainX = subX(:,:,:,cv.training);
            testX = subX(:,:,:,cv.test);
            trainY = Y(cv.training);
            testY = Y(cv.test);

            %{
                maxPooling2dLayer([3 3],'Stride',2)
                convolution2dLayer([3 3],3)
                reluLayer 
            %}
            %% 11 - 5
            layers = [ ...
                imageInputLayer(szX(1:2),'Normalization','none')
                convolution2dLayer([18 18],5)
                reluLayer
                
                fullyConnectedLayer(1)
                regressionLayer];
            vw = 'training-progress';
            sc = 'every-epoch';
            options = trainingOptions('sgdm', ...
                'L2Regularization',.001,...
                'LearnRateSchedule','piecewise',...
                'LearnRateDropPeriod',50,...
                'LearnRateDropFactor',1,...
                'MiniBatchSize',32,...
                'Shuffle',sc,...
                'InitialLearnRate',.000001,...
                'ValidationData',{testX,testY'},...
                'MaxEpochs',1000, ...
                'ValidationPatience',10,...
                'ExecutionEnvironment','parallel',...
                'Plots',vw);

            [countNet] = trainNetwork(trainX,trainY',layers,options);
            
%% brute force count
convolutionSZ = 10:5:25;
filterNumber = 3:2:21;
traininfo = {};
LOOP = 3;
toPlot = [];
close all force
figure;
for f = 1:numel(filterNumber)
    for c = 1:numel(convolutionSZ)
        loopValue = [];
        for loop = 1:LOOP
           
            subX = [];
            trainY = Y;

            for e = 1:size(X,3)
                subX(:,:,e) = imresize(X(:,:,e),.25);
            end

            
            szX = size(subX);
            subX = reshape(subX,[szX(1) szX(2) 1 szX(3)]);

            cv = cvpartition(size(subX,4),'HoldOut',.2);
            trainX = subX(:,:,:,cv.training);
            testX = subX(:,:,:,cv.test);
            trainY = Y(cv.training);
            testY = Y(cv.test);

            layers = [ ...
                imageInputLayer(szX(1:2),'Normalization','none')
                convolution2dLayer([convolutionSZ(c) convolutionSZ(c)],filterNumber(f))
                reluLayer

                fullyConnectedLayer(1)
                regressionLayer];

            %{
            layers = [ ...
                imageInputLayer(szX(1:2),'Normalization','none')
                convolution2dLayer([11 11],13)
                reluLayer

                maxPooling2dLayer([3 3],'Stride',2)
                convolution2dLayer([3 3],3)
                reluLayer 

                fullyConnectedLayer(1)
                regressionLayer];
            %}

            [trainY_z,mu,sigma] = zscore(trainY);
            %{
            options = trainingOptions('sgdm', ...
                'MiniBatchSize',32,...
                'Shuffle','every-epoch',...
                'InitialLearnRate',.001,...
                'ValidationData',{xTest,ZyTest},...
                'MaxEpochs',1000, ...
                'ExecutionEnvironment','parallel',...
                'Plots','training-progress');
            %}

            vw = 'training-progress';
            %vw = 'none';
            
            %sc = 'every-epoch';
            sc = 'once';
            options = trainingOptions('sgdm', ...
                'MiniBatchSize',32,...
                'Shuffle',sc,...
                'InitialLearnRate',.0001,...
                'ValidationData',{testX,testY'},...
                'MaxEpochs',1000, ...
                'ValidationPatience',3,...
                'ExecutionEnvironment','parallel',...
                'Plots',vw);

            [countNet,traininfo{c,f,loop}] = trainNetwork(trainX,trainY',layers,options);
            
            kidx = ~isnan(traininfo{c,f,loop}.ValidationRMSE);
            loopValue(loop) = min(traininfo{c,f,loop}.ValidationRMSE(kidx));
        end
        close all force
        toPlot(c,f) = min(loopValue);
        plot(toPlot);
        drawnow
    end
end
%% test the count net
close all force
squareSZ = [1000 1000];
for e = 2:10
    videoFile = FileList{e};
    skipT = 10;
    runningCount = [];
    cnt = 1;
    for t = 1:skipT:1500
        t
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % read the image(s)
        I = readSlice(videoFile,experimentFrame(e) + t,0);
        % get the petri mask from first pass methed(s)
        petriMask = logical(maskStack(:,:,e));
        % get the bounding box for the petri plate
        pR = regionprops(petriMask);
        % crop mask and image and background
        I = imcrop(I,pR.BoundingBox);
        M = imcrop(petriMask,pR.BoundingBox);
        BK = imcrop(background(:,:,e),pR.BoundingBox);
        I = (I - BK).*M;
        %I = I.*M;
        %Io = I;
        Io = imresize(I,squareSZ);
        
        I = imresize(Io,per);
        tQ = im2colF(I,[smSZ smSZ],[1 1]);
        tQ = bsxfun(@times,flt(:),tQ);
        tQ = sort(tQ,1);


        tmp = imresize(double(edge(Io)),per);
        tQ2 = im2colF(tmp,[smSZ smSZ],[1 1]);
        tQ2 = bsxfun(@times,flt(:),tQ2);
        tQ2 = sort(tQ2,1);
        
        I = PCA_REPROJ_T([tQ;tQ2],E,U);
        I = sim(cNet,I);
        I = col2im(I,[21 21],size(tmp));
        
        
        %runningCount(cnt) =  (sS*countNet.predict(I))+sU;
        runningCount(cnt) =  countNet.predict(I);
        cnt = cnt + 1;
        
        plot(runningCount,'k');hold on
        plot(mean(runningCount)*ones(size(runningCount)),'g--');
        plot(median(runningCount)*ones(size(runningCount)),'b--');
        plot(totalCounts(e)*ones(size(runningCount)),'k--');
        axis([1 numel(runningCount)+1 0 50]);
        hold off
        drawnow
        
    end
    plot(runningCount,'k');hold on
    plot(mean(runningCount)*ones(size(runningCount)),'g--');
    plot(median(runningCount)*ones(size(runningCount)),'b--');
    plot(totalCounts(e)*ones(size(runningCount)),'k--');
    axis([1 numel(runningCount)+1 0 50]);
    hold off
    drawnow
    waitforbuttonpress
end
%% loop plot pretty each of the datasets
X = [];
Y = [];
for e = 1:numel(NUM_BEST)


    
    tmC = timeCounts(e,:);

    ON_BEST = NUM_BEST{e};
    ON = NUM_IN{e};
    VN = NUM_VEL{e};

    t = numel(ON);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tvec = (1:t) - 1;
    %plot(tvec,ON,'k');hold on
    %plot(tvec,VN,'r');
    %plot(tvec,ON_BEST,'c');



    av = .5*mean(ON_BEST)+.5*mean(ON);




    LEG{1} = 'Manual Total Count';
    LEG{2} = 'Automated Total Count';
    LEG{3} = 'Manual Paralyzed Count ';
    LEG{4} = 'Automated Non-Moving Count';

    plot(tvec,totalCounts(e)*ones(size(ON)),'k','LineWidth',1);

    hold on
    %plot(tvec,ON_BEST)
    %plot(tvec,mean(ON_BEST)*ones(size(ON)),'c--')
    %plot(tvec,mean(ON)*ones(size(ON)),'k--')
    plot(tvec,(av)*ones(size(ON)),'r--');


    plot(timeHeader*10,tmC,'k','LineWidth',1);

    smoothAmount = 11;
    notMoving1 = .5*(ON_BEST+ON) - VN;
    %plot(tvec,notMoving1,'b-','LineWidth',.5)
    smoothNotMoving1 = imfilter(notMoving1,ones(1,smoothAmount)/smoothAmount,'replicate');
    %plot(tvec,smoothNotMoving1,'b','LineWidth',1);


    notMoving2 = ON_BEST - VN;
    %plot(tvec,notMoving2,'c-','LineWidth',.5)
    smoothNotMoving2 = imfilter(notMoving2,ones(1,smoothAmount)/smoothAmount,'replicate');
    plot(tvec,smoothNotMoving2,'r','LineWidth',1);

    %{
    notMoving3 = mean(ON_BEST)*ones(size(ON_BEST)) - VN;
    %plot(tvec,notMoving2,'c-','LineWidth',.5)
    smoothNotMoving3 = imfilter(notMoving3,ones(1,smoothAmount)/smoothAmount,'replicate');
    plot(tvec,smoothNotMoving2,'c','LineWidth',1);
    %}


    smoothAmount = 5;
    smoothV = imfilter(VN,ones(1,smoothAmount)/smoothAmount,'replicate');
    notMoving3 = mean(ON_BEST)*ones(size(ON_BEST)) - smoothV;
    %plot(tvec,notMoving2,'c-','LineWidth',.5)
    smoothNotMoving3 = imfilter(notMoving3,ones(1,smoothAmount)/smoothAmount,'replicate');
    %plot(tvec,smoothNotMoving3,'m','LineWidth',1);



    sam = 1:50:numel(tvec);
    %sam = timeHeader*10;
    %plot(tvec(sam),smoothNotMoving2(sam),'r')
    %plot(tvec(sam),smoothNotMoving2(sam),'ro')

   % plot(timeHeader*10,tmC,'k','LineWidth',1);
    %plot(timeHeader*10,tmC,'k*');

    hold off
    axis([0 t -10 50]);
    xlabel(['Time (frames)']);
    ylabel(['Counts']);
    legend(LEG)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    waitforbuttonpress
    pause(.5)

    tX = [smoothNotMoving3(sam)',smoothNotMoving2(sam)',smoothNotMoving1(sam)'];
    tY = [tmC'];


    X = [X;tX];
    Y = [Y;tY];
   
end
                
%% look at vStack and iStack
data = [velStack(:) iStack(:)];
ridx = randperm(size(data,1));
%%
sam = 100000000;
sub = data(ridx(1:sam),:);
%% remove high noise
sub(sub(:,1) > 10,:) = [];
%% 2D scatter
plot(sub(:,1),sub(:,2),'.')
%% ksdensity
close all
[f,xi] = ksdensity(sub(:,1));
f = f / sum(f);
plot(xi,log(f))
xi = linspace(xi(1),xi(end),1000);
xi = linspace(0,10,10000);
[f,xi] = ksdensity(sub(:,1),xi);

plot(xi,log(f));
%%
f(1) = [];
xi(1) = [];
f(end) = [];
xi(end) = [];
plot(xi,log(f));
%f = f / sum(f);
%%
close all
clear p
d1 = @(x,p)pdf('Exponential',x,p(1));
d2 = @(x,p)pdf('Normal',x,p(1),p(2));
n1 = @(x,p)d1(x,p)/sum(d1(x,p));
n2 = @(x,p)d2(x,p)/sum(d2(x,p));
g1 = @(x,p)p(1)*n1(x,p(2:end));
g2 = @(x,p)p(1)*n2(x,p(2:end));
dt = @(x,p)g1(x,p(1:2)) + g2(x,p(3:5));
%dt = @(x,p)n1(x,p(2));
contrast = @(p)norm(log(dt(xi,p)) - log(f));
%contrast = @(p)-log(sum(normpdf(log(dt(xi,p)) - log(f),0,1)));
ip = [.9 .18 .1 4 100];
y = dt(xi,ip);
%%
MX = 2;
data = sub(1:1000000,1);
data(data > MX) = [];

xi = linspace(0,MX,1000);
[f,xi] = ksdensity(data,xi);

clear x
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
c = @(x,p)exp((p(1)/p(2))*(1-x.^p(2))).*x.^(p(2)-2).*(p(1)*x.^p(2) - p(2) + 1);
lb = [0.00001;0];
ub = [10000;.9999];
ip = [3 1];
mx = min(data(:));
t = @(x)(x-mx+1);
%plot(xi,c(xi,[3 1]))
%waitforbuttonpress
%}
pd = fitdist(data(:),'Weibull');
pd = fitdist(data(:),'Exponential');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = @(x,p)(p(2)/p(1))*((x*p(1)^-1).^(p(2)-1)).*(exp(-(x*p(1)^-1).^p(2)));
lb1 = [0;0];
ub1 = [10000;100000];
%ip1 = [pd.A pd.B];
ip1 = [1 1];
t = @(x)x;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q = @(x,p)normpdf(x,p(1),p(2));
lb2 = [0;1];
ub2 = [10;100];
ip2 = [7 2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = @(x,p)p(1)*(1-exp(-p(2)*x));
t = @(x)(x - x(1));
m = @(p)norm(n(xi,p) - t(-log(f)));
m = @(p)-sum(log(normpdf((n(xi,p) - t(-log(f))),0,1)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mw = @(x,p)((2/pi)^.5)*((x.^2).*exp((-x.^2).*(2*p(1)^2)^-1))*p(1)^-3;
mw = @(x,p)(x*p(1)^-2).*(exp((-x.^2)*(2*p(1)^2).^-1));
% dagum
mw = @(x,p)(p(1)*p(2)*(x.^-1)).*(((x*p(3)^-1).^(p(1)*p(2))).*((x*p(3).^-1).^(p(1)) + 1).^-(p(2)+1));
%mw = @(x,p)(p(1)*p(2))*((x.^(p(1)-1)).*((1+x.^p(1)).^(p(2)+1)).^-1)
%mw = @(x,p)((x.^(p(1)-1)).*((1+x).^(-p(1)-p(2))))*beta(p(1),p(2))^-1;
m = @(p)-sum(log(normpdf((mw(xi+.0000000000001,p) - f))));
ip = .001;
ip = .01;
% dagum
ip = [1 2 .005];
lb = [0;0;0];
ub = [1000,1000,1000];


%ip = [1.1 .2];
%lb = [0;0];
%ub = [1000,1000];

%k = @(x,p)p(1)*c(x,2:3) + p(4)*q(x,p(5:end));
%lb = [0;lb1;0;lb2];
%ub = [1;ub1;1;ub2];
%ip = [.99 ip1 .01 ip2];

%pd = fitdist(sub(:,1),'Weibull');
%ip = [pd.A pd.B];


contrast = @(p)-sum(log(mw(t(data),p)));
%contrast = @(p)-sum(log(c(t(data),p)));
%contrast = @(p)-sum(log(k(data,p)));


%contrast = @(p)-sum(log(c(xi,p)));
%contrast(ip)
%contrast(ip*1.001)
%contrast(ip*.99)

ops = optimoptions('patternsearch','Display','iter','UseParallel',true,...
                    'MaxIterations',1000000,'MeshTolerance',10^-10,'StepTolerance',10^-10);
Aeq = [1 0 0 1 0 0];
beq = 1;

Aeq = [];
beq = [];
w = patternsearch(m,ip,[],[],Aeq,beq,lb,ub,[],ops);
%plot(xi,log(pdf(pd,xi)),'g-');hold on
%plot(xi,c(xi,[pd.A pd.B]),'k--')
%y = c(xi,w);
y = mw(xi,w);
plot(xi,log(y),'b-');hold on
plot(xi,log(f),'r--');
plot(xi,log(f-y),'k')
hold on
%plot(xi,t(-log(f)),'r');
%u = fminsearch(m,[-.5,mean(t(xi((end-100):end)))]);
%plot(xi,n(xi,u),'b');
%figure;
%plot(xi,-log(f));hold on
%plot(xi,n(xi,u)-log(f(1)))
%figure;
%plot(xi,exp(log(f)),'r');hold on
%plot(xi,exp(-(n(xi,u)-log(f(1)))),'k')
%%
figure;

MX = .1;
data = sub(1:100000,1);
data(data > MX) = [];

pd = fitdist(data,'Weibull');

xi = linspace(0,MX,1000);
[f,xi] = ksdensity(data,xi);
plot(xi,f,'r');hold on
plot(xi,pdf(pd,xi),'k')
%%
close all
%pd = fitdist(sub(:,1),'Exponential');
pd = fitdist(sub(:,1),'Weibull');

plot(xi,-log(pdf(pd,xi)))
hold on
plot(xi,f,'r')
%%
ops = optimoptions('patternsearch','Display','iter','UseParallel',true,...
                    'MaxIterations',1000000,'MeshTolerance',10^-10,'StepTolerance',10^-10);
Aeq = [1 0 1 0 0];
beq = 1;


lb = zeros(numel(ip),1);
lb(2) = -10;
lb(4) = 3;
ub = [1 1000 .1 10 1000];
w = patternsearch(contrast,ip,[],[],Aeq,beq,lb,ub,[],ops);
sy1 = g1(xi,w(1:2));
sy1 = sy1 / sum(sy1);

y = dt(xi,w);
y = y / sum(y);
t = @(x)log(x);
t = @(x)x;
close all
plot(xi,t(sy1),'k');hold on
plot(xi,t(f),'r')
plot(xi,t(y),'b');hold on
%%
y = dt(xi,w);
plot(xi,log(y),'k')
hold on
plot(xi,log(f),'r')
%%
plot(xi,log(f))
f = log(f);
f = bindVec(f);
f = f / sum(f);
plot(xi,log(f));
%% 
close all
for e = 1:numel(V)
    sig = V{e};
    sig = imfilter(sig,fspecial('average',11),'replicate');
    plot(sig)
    hold all
    axis([0 500 0 10])
end
%% read frame(s) slice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
videoFile = FileList{1};
M = VideoReader(videoFile);
T = M.NumberofFrames;
H = M.Height;
W = M.Width;
[slice] = readSlice(videoFile,100,10);
per = .25;
[nslice] = resizeSlice(slice,per);
%% construct the tracker
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I(:,:,1)      is the first image
% I(:,:,2)      is the second image
% para.P        is the point to track
% para.domain   is the tracking domain
% para.RADIUS   is the clipping window
% para.init_T   is the init vector for transformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
movieSubSelect = 5;
% get random slice from movie
[trackSlice,trackMovieSelect] = randomSlice(FileList(movieSubSelect),5,1);
% get the mask(s)
trackMask = generateAssociatedMaskSlice(movieSubSelect,maskStack);
% generate select points
[trackPointMask] = generateSelectPoints(trackMask,[10 10]);
[trackPoints(:,2),trackPoints(:,1)] = find(trackPointMask);
% get the temp track slice
toTrack1 = squeeze(trackSlice(:,:,:,1));
toTrack2 = squeeze(trackSlice(:,:,:,2));
% make the domain
RAD = 45;NP = [RAD 100];
[rad,theta] = ndgrid(linspace(0,RAD,NP(1)),linspace(-pi,pi,NP(2)));
trackD = [rad(:).*cos(theta(:)),rad(:).*cos(theta(:))];
% create function to generate the transformation at P
fT = @(p,globData)[[eye(2) p'];[0 0 1]];

P = trackPoints(1,:);

mov = squeeze(trackSlice);
close all
for e = 1:size(mov,3)
    imshow(mov(:,:,e),[]);
    drawnow
end
[c,r,V] = impixel(mov(:,:,1));

%%
flowTrack = [c r];
P = [c r];
dP = 20;
dN = 5;
[p1,p2] = ndgrid(linspace(P(1)-dP,P(1)+dP,dN),linspace(P(2)-dP,P(2)+dP,dN));
P = [p1(:),p2(:)];
clear trackStack
for p = 1:size(P,1)

    initP1 = funcQ(P(p,:),fT,[]);
    initP1.getT([]);

    trackStack(p) = initP1;
end

restrictionSize = 51;
RAD = 25;
[y,x] = ndgrid(linspace(-RAD,RAD,2*RAD),linspace(-RAD,RAD,2*RAD));
domain.sz = size(x);
domain.domain = [x(:),y(:),ones(size(x(:)))];


RAD = 35;NP = [RAD 75];
[rad,theta] = ndgrid(linspace(0,RAD,NP(1)),linspace(-pi,pi,NP(2)));
domain.sz = size(rad);
domain.domain = [rad(:).*cos(theta(:)),rad(:).*sin(theta(:)),ones(size(rad(:)))];


opticFlow1 = opticalFlowFarneback('NeighborhoodSize',15);

for t = 1:(size(trackSlice,4)-1)

    toTrack1 = squeeze(trackSlice(:,:,:,t));
    toTrack2 = squeeze(trackSlice(:,:,:,t+1));
    
    
    stackSlice = trackStack(t,:);
    nextSlice = funcQ(rand([1,2]),fT,[]);
    nextSlice = repmat(nextSlice,size(stackSlice));


    flowField = estimateFlow(opticFlow1,toTrack1);


    
    dX = [ba_interp2(flowField.Vx,flowTrack(end,1),flowTrack(end,2)),
    ba_interp2(flowField.Vy,flowTrack(end,1),flowTrack(end,2))];

    flowField = [flowField;flowField(end,:) + dX];


    for pr = 1:size(trackStack,2)
        pr;
        tmpP = stackSlice(pr).P;
        [subI1,newP,upperLeft] = cropAtP(toTrack1,tmpP,restrictionSize);
        stackSlice(pr).globData = subI1;
        stackSlice(pr).setP(newP);
        
        [subI2,newP,upperLeft] = cropAtP(toTrack2,tmpP,restrictionSize);
        nextPStatic = funcQ(newP,fT,subI2);
        nextPStatic.resetT();
        
        nextP = stackSlice(pr).morph(nextPStatic,domain);
        nextP.P = nextP.P + upperLeft;
        nextP.getT([]);
        
        nextSlice(pr) = nextP;
        
    end

    trackStack = [trackStack;nextSlice];
     

    imshow(toTrack2,[]);
    hold on;

    for pr = 1:size(trackStack,2)
        plot(trackStack(t+1,pr).P(1),trackStack(t+1,pr).P(2),'r.');
    end
    drawnow


    hold off
    drawnow

end
%%
close all
for t = 1:size(trackStack,1)
    toTrack2 = squeeze(trackSlice(:,:,:,t+1));
    imshow(toTrack2,[]);
    hold on
    for pr = 1:size(trackStack,2)
        plot(trackStack(t+1,pr).P(1),trackStack(t+1,pr).P(2),'r.');
    end
    hold off
    drawnow
waitforbuttonpress
end
%%


%{
%%
% generate domain for tracking
RAD = 45;NP = [RAD 100];
[rad,theta] = ndgrid(linspace(0,RAD,NP(1)),linspace(-pi,pi,NP(2)));
trackD = [rad(:).*cos(theta(:)),rad(:).*cos(theta(:))];
init_T = [[1 0 0];[0 1 0]];
trackPara.P = trackPoints(1,:);
trackPara.domain = trackD;
trackPara.RADIUS = RAD+20;
trackPara.init_T = init_T(:)';
trackPara.externalInterp = true;
trackPara.threshold = .001;
toTrack = squeeze(trackSlice(:,:,:,1:2));
[T] = miniTrack(toTrack,trackPara);
%}
%% build up pca space first for sorted and resized
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get 100 random slices
[rslice,selectedTrainIDX] = randomSlice(FileList,0,100);
% get the mask(s)
mslice = generateAssociatedMaskSlice(selectedTrainIDX,maskStack);
% cat the masks
rslice = cat(3,rslice,mslice);
% resize the slices and masks
[rslice] = resizeSlice(rslice,per);

szV = [31 31];
skV = [3 3];
func = @(x)sort(x,1);
[v] = vectorizeSlices(rslice,szV,skV,func);
szV = size(v);
v = reshape(v,[szV(1) prod(szV(2:end))]);

%%
[srtU,srtE,srtL] = PCA_FIT_FULL_Tws(v,3);
%% build up the pca space with sampler
[rslice,selectedTrainIDX] = randomSlice(FileList,0,100);
%% build up w-slices
Nsam = 100;
[rslice,selectedTrainIDX] = randomSlice(FileList,2,Nsam);
[rslice] = resizeSlice(rslice,per);
szV = [31 31];
skV = [1 1];
func = @(x)sort(x,1);
func2 = @(x)PCA_REPROJ_T(sort(x,1),srtE,srtU,0);
[v] = vectorizeSlices(rslice,szV,skV,func2);
v = permute(v,[1 3 2 4]);
szV = size(v);
v = reshape(v,[prod(szV(1:2)) szV(3:4)]);
%% try CNN - step 1 sample slices
%%%%%%%%%%%%%%%%%%%%%%%%
c = cvpartition(numel(FileList),'HoldOut',.2);
Nsam = 100;
per = .35;
WIN = 5;
trainIdx = find(c.training);
testIdx = find(c.test);
[rslice,selectedTrainIDX,FR] = randomSlice(FileList(trainIdx),WIN,Nsam);
[rslice] = resizeSlice(rslice,per);
rmidx = find(FR < (experimentFrame(trainIdx(selectedTrainIDX)) + WIN));
rslice(:,:,:,:,rmidx) = [];
selectedTrainIDX(rmidx) = [];
%% gather test set
Nsam = 40;
[testSlice,selectedTestIDX,FR] = randomSlice(FileList(testIdx),WIN,Nsam);
[testSlice] = resizeSlice(testSlice,per);
rmidx = find(FR < (experimentFrame(testIdx(selectedTestIDX)) + WIN));
testSlice(:,:,:,:,rmidx) = [];
selectedTestIDX(rmidx) = [];
%% gather mask data - training
%%%%%%%%%%%%%%%%%%%%%%%%
close all
mslice = generateAssociatedMaskSlice(trainIdx(selectedTrainIDX),maskStack);
[mslice] = resizeSlice(mslice,per);
for e = 1:size(mslice,5)
    tmp = mslice(:,:,1,1,e) > .5;
    velR = regionprops(tmp,'BoundingBox');
    bBOX(e,:) = velR.BoundingBox;
end
boxSZ = max(bBOX(:,3:4),[],1);
TRAINbslice = [];
% for each sample
for e = 1:size(mslice,5)
    tmp = mslice(:,:,1,1,e) > .5;
    velR = regionprops(tmp,'BoundingBox');
   	box = velR.BoundingBox;
    box(3:4) = boxSZ;
    tmpM = imcrop(squeeze(mslice(:,:,1,1,e)),box);
    for sl = 1:size(rslice,4)
        tmp = imcrop(squeeze(rslice(:,:,1,sl,e)),box);
        TRAINbslice(:,:,1,sl,e) = tmp.*tmpM;
    end
    imshow(squeeze(TRAINbslice(:,:,1,1,e)),[]);
    drawnow
end
%% gather mask data - test
%%%%%%%%%%%%%%%%%%%%%%%%
close all
mslice = generateAssociatedMaskSlice(testIdx(selectedTestIDX),maskStack);
[mslice] = resizeSlice(mslice,per);
for e = 1:size(mslice,5)
    tmp = mslice(:,:,1,1,e) > .5;
    velR = regionprops(tmp,'BoundingBox');
    bBOX(e,:) = velR.BoundingBox;
end
boxSZ = max(bBOX(:,3:4),[],1);
TESTbslice = [];
% for each sample
for e = 1:size(mslice,5)
    tmp = mslice(:,:,1,1,e) > .5;
    velR = regionprops(tmp,'BoundingBox');
   	box = velR.BoundingBox;
    box(3:4) = boxSZ;
    tmpM = imcrop(squeeze(mslice(:,:,1,1,e)),box);
    for sl = 1:size(rslice,4)
        tmp = imcrop(squeeze(rslice(:,:,1,sl,e)),box);
        TESTbslice(:,:,1,sl,e) = tmp.*tmpM;
    end
    imshow(squeeze(TESTbslice(:,:,1,1,e)),[]);
    drawnow
end
%%
Xplay = squeeze(TRAINbslice);
szP1 = size(Xplay);
Xplay = permute(Xplay,[3 1 2 4]);
szP2 = size(Xplay);
Xplay = reshape(Xplay,[szP2(1) prod(szP2(2:end))]);
Sp = PCA_REPROJ_T(Xplay,Ep,Up);
szP3 = szP2;
szP3(1) = cn;
Sp = reshape(Sp,szP3);
Sp = ipermute(Sp,[3 1 2 4]);
xTrain = Sp;
%%
Xplay = squeeze(TESTbslice);
szP1 = size(Xplay);
Xplay = permute(Xplay,[3 1 2 4]);
szP2 = size(Xplay);
Xplay = reshape(Xplay,[szP2(1) prod(szP2(2:end))]);
Sp = PCA_REPROJ_T(Xplay,Ep,Up);
szP3 = szP2;
szP3(1) = cn;
Sp = reshape(Sp,szP3);
Sp = ipermute(Sp,[3 1 2 4]);
xTest = Sp;
%% gather
%%%%%%%%%%%%%%%%%%%%%%%
for e = 1:size(TRAINbslice,5)
    tmp = TRAINbslice(:,:,:,:,e);
    szT = size(tmp);
    tmp = reshape(tmp,[prod(szT(1:2)) szT(3:end)]);
    offset = mean(tmp,1);
    tmp = bsxfun(@minus,tmp,offset);
    tmp = reshape(tmp,szT);
    TRAINbslice(:,:,:,:,e) = tmp;
end

for e = 1:size(TESTbslice,5)
    tmp = TESTbslice(:,:,:,:,e);
    szT = size(tmp);
    tmp = reshape(tmp,[prod(szT(1:2)) szT(3:end)]);
    offset = mean(tmp,1);
    tmp = bsxfun(@minus,tmp,offset);
    tmp = reshape(tmp,szT);
    TESTbslice(:,:,:,:,e) = tmp;
end
%%
% prepare the train data
%%%%%%%%%%%%%%%%%%%%%
uTrain = squeeze(mean(TRAINbslice,4));
uTrain = squeeze(max(TRAINbslice,[],4));
szU = size(uTrain);
szU = [szU(1:2) 1 szU(3)];
uTrain = reshape(uTrain,szU);
sTrain = squeeze(std(TRAINbslice,1,4));
sTrain = reshape(sTrain,szU);
zTrain = zeros(size(sTrain));
zTrain = mean(abs(diff(TRAINbslice,1,4)),4);
zTrain = reshape(zTrain,szU);
xTrain = cat(3,uTrain,sTrain,zTrain);
for e = 1:size(xTrain,4)
    for k = 1:size(xTrain,3)
        xTrain(:,:,k,e) = bindVec(xTrain(:,:,k,e));
    end
end

% prepare test data
%%%%%%%%%%%%%%%%%%%%%
uTest = squeeze(mean(TESTbslice,4));
uTest = squeeze(max(TESTbslice,[],4));
szU = size(uTest);
szU = [szU(1:2) 1 szU(3)];
uTest = reshape(uTest,szU);
sTest = squeeze(std(TESTbslice,1,4));
sTest = reshape(sTest,szU);
zTest = zeros(size(sTest));
zTest = mean(abs(diff(TESTbslice,1,4)),4);
zTest = reshape(zTest,szU);
xTest = cat(3,uTest,sTest,zTest);
for e = 1:size(xTest,4)
    for k = 1:size(xTest,3)
        xTest(:,:,k,e) = bindVec(xTest(:,:,k,e));
    end
end
%%

yTrain = totalCounts(trainIdx(selectedTrainIDX));
yTest = totalCounts(testIdx(selectedTestIDX));
szX = size(xTrain);

layers = [ ...
    imageInputLayer(szX(1:3),'Normalization','none')
    convolution2dLayer([21 21],5)
    reluLayer
    fullyConnectedLayer(1)
    regressionLayer];


%{
layers = [ ...
    imageInputLayer(szX(1:3),'Normalization','none')

    convolution2dLayer([21 21],5)
    reluLayer


    maxPooling2dLayer(2,'Stride',2)

    convolution2dLayer(5,5,'Padding',1)
    reluLayer

    fullyConnectedLayer(1)
    regressionLayer];
%}
    
    
    
%{
layers = [ ...
    imageInputLayer(szX(1:3),'Normalization','none')
    convolution2dLayer([31 31],15)
    reluLayer
    fullyConnectedLayer(1)
    regressionLayer];
%}

layers = [ ...
    imageInputLayer(szX(1:3),'Normalization','none')
    convolution2dLayer([15 15],3,'Stride',7)
    reluLayer
    fullyConnectedLayer(1)
    regressionLayer];



[ZyTrain,mu,sigma] = zscore(yTrain);
ZyTest = bsxfun(@times,bsxfun(@minus,yTest,mu),sigma.^-1);
%ZyTest = yTest;
%ZyTrain = yTrain;
options = trainingOptions('sgdm', ...
    'MiniBatchSize',32,...
    'Shuffle','every-epoch',...
    'InitialLearnRate',.001,...
    'ValidationData',{xTest,ZyTest},...
    'MaxEpochs',1000, ...
    'ExecutionEnvironment','parallel',...
    'Plots','training-progress');

countNet = trainNetwork(xTrain,ZyTrain,layers,options);
%% test count net
Nsam = 1;
[testSlice,testIDX2] = randomSlice(FileList(testIdx),2,Nsam);
testMaskSlice = generateAssociatedMaskSlice(testIDX2,maskStack);
testSlice = resizeSlice(testSlice,per);
[testMaskSlice] = resizeSlice(testMaskSlice,per);
countSlice(countNet,testSlice,testMaskSlice,boxSZ)
totalCounts(testIdx(testIDX2));
%% generate short read basis vectors
close all
Xplay = squeeze(TRAINbslice);
szP1 = size(Xplay);
Xplay = permute(Xplay,[3 1 2 4]);
szP2 = size(Xplay);
Xplay = reshape(Xplay,[szP2(1) prod(szP2(2:end))]);
toUse = find(~all(Xplay==0,1));
cn = 4;

mini = min(Xplay,[],1);
nXplay = bsxfun(@minus,Xplay,mini);
maxi = max(nXplay,[],1);
nXplay = bsxfun(@times,nXplay,maxi.^-1);



[Up,Ep,Lp]  = PCA_FIT_FULL_Tws(Xplay(:,toUse),cn);
[Upn,Epn,Lpn]  = PCA_FIT_FULL_Tws(nXplay(:,toUse),cn);

Sp = PCA_REPROJ_T(Xplay,Ep,Up);
szP3 = szP2;
szP3(1) = cn;
Sp = reshape(Sp,szP3);
Sp = ipermute(Sp,[3 1 2 4]);

%% generate small read basis vectors
SAM_SP = [];
for movSEL = 1:size(TRAINbslice,5)

movSEL

    close all
    %movSEL = 3;
    gz = 5;
    tmp = squeeze(TRAINbslice(:,:,:,1,movSEL));
    scanSliceMask = (tmp(:,:,1) ~= 0);
    scanSliceMask = imerode(scanSliceMask,strel('disk',31,0));
    [g1,g2] = gradient(tmp);
    [m1,m2] = find(scanSliceMask);
    midx = find(scanSliceMask);
    [n1,n2] = ndgrid(-gz:gz,-gz:gz);
    n = [n1(:) n2(:) ones(size(n1(:)))];
    imshow(tmp,[]);
    hold on
    G1 = ba_interp2(g1,m2,m1);
    G2 = ba_interp2(g2,m2,m1);
    TAN = [G1 G2];
    TAN = bsxfun(@times,TAN,sum(TAN.*TAN,2).^-.5);
    NOR = [TAN(:,2) -TAN(:,1)];
    imshow(tmp,[]);
    hold on
    quiver(m2,m1,TAN(:,1),TAN(:,2));
    MID = [m2 m1];
    close all
    affine = [];
    for e = 1:numel(m1)
        %imshow(tmp,[]);
        %hold on
        affine(:,:,e) = [[TAN(e,:)' NOR(e,:)' MID(e,:)'];[0 0 1]];
        %nt = (affine(:,:,e)*n')';
        %plot(nt(:,1),nt(:,2),'.')
        %hold off
        %drawnow
    end

    affine = permute(affine,[3 2 1]);
    nt = mtimesx(affine,n,'T');

    nt = nt(:,:,1:2);
    nt = permute(nt,[2 1 3]);
    szN = size(nt);

    nt = reshape(nt,[prod(szN(1:2)) szN(3)]);
    close all
    imshow(tmp,[]);
    hold on
    plot(nt(1:121,1),nt(1:121,2),'.')
    hold off
    F = ba_interp2(tmp,nt(:,1),nt(:,2));
    F = reshape(F,[szN(1) szN(2)]);


    Fc = PCA_REPROJ_T(F,Eq,Uq,0);
    
    Z = zeros([size(tmp) cn2]);
    for k = 1:size(Z,3)
        tmpK = Z(:,:,k);
        tmpK(midx) = Fc(k,:);
        Z(:,:,k) = tmpK;
    end

    SAM_SP = [SAM_SP F];

end
cn2 = 3;
[Uq,Eq,Lq]  = PCA_FIT_FULL_Tws(SAM_SP,cn2);
%% scan movie
masterSpace = [];
for movieSEL = 1:numel(FileList)
    %movieSEL = 4;
    frameSKIP = 150;
    videoFile = FileList{movieSEL};
    M = VideoReader(videoFile);
    T = M.NumberofFrames;
    TOT = T - WIN;

    scanSliceMask = generateAssociatedMaskSlice(movieSEL,maskStack);
    [scanSliceMask] = resizeSlice(scanSliceMask,per);
    maskIDX = find(scanSliceMask ~= 0);

    for fr = experimentFrame(movieSEL):frameSKIP:TOT
        [scanSlice,szComplex] = readSlice(videoFile,fr,WIN);
        [scanSlice] = resizeSlice(scanSlice,per);
        scanSlice = bsxfun(@times,scanSlice,scanSliceMask);



        centerImage = squeeze(scanSlice(:,:,1,(end-1)/2));
        [spaceChunk,spaceMaskIdx] = spaceBLOCK(centerImage,5,Eq,Uq,cn2);

        spaceChunk_vec = spaceChunk;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % color normalize
        for k = 1:size(spaceChunk,3)
            tmp = spaceChunk(:,:,k);
            tmp(spaceMaskIdx) = bindVec(tmp(spaceMaskIdx));
            %tmp(spaceMaskIdx) = tmp(spaceMaskIdx) - mean(tmp(spaceMaskIdx));
            spaceChunk(:,:,k) = tmp;
        end
        %imshow(spaceChunk,[]);
        %drawnow
        %waitforbuttonpress


        Xplay = squeeze(scanSlice);
        szP1 = size(Xplay);
        Xplay = permute(Xplay,[3 1 2 4]);
        szP2 = size(Xplay);
        Xplay = reshape(Xplay,[szP2(1) prod(szP2(2:end))]);

        szP3 = szP2;
        szP3(1) = cn-1;

        Sp = zeros(szP3);
        Sp(:,maskIDX) = PCA_REPROJ_T(Xplay(:,maskIDX),Ep(:,2:end),Up);

        Sp = reshape(Sp,szP3);
        Sp = ipermute(Sp,[3 1 2 4]);

        Spn = zeros(szP3);

        mini = min(Xplay,[],1);
        nXplay = bsxfun(@minus,Xplay,mini);
        maxi = max(nXplay,[],1);
        nXplay = bsxfun(@times,nXplay,maxi.^-1);
        Spn(:,maskIDX) = PCA_REPROJ_T(nXplay(:,maskIDX),Epn(:,2:end),Upn);
        Spn = reshape(Spn,szP3);
        Spn = ipermute(Spn,[3 1 2 4]);




        timeChunk = Sp;
        timeChunk_vec = timeChunk;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % color normalize
        for k = 1:size(timeChunk,3)
            tmp = timeChunk(:,:,k);
            tmp(maskIDX) = bindVec(tmp(maskIDX));
            %tmp(maskIDX) = tmp(maskIDX) - mean(tmp(maskIDX));
            timeChunk(:,:,k) = tmp;
        end

        timeChunk_n = Spn;
    %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % color normalize
        for k = 1:size(timeChunk_n,3)
            tmp = timeChunk_n(:,:,k);
            tmp(maskIDX) = bindVec(tmp(maskIDX));
            tmp(maskIDX) = tmp(maskIDX) - mean(tmp(maskIDX));
            timeChunk_n(:,:,k) = tmp;
        end
    %}
        %imshow(timeChunk,[]);
        %drawnow
        %waitforbuttonpress

        square1 = ~all(centerImage == 0,1);
        square2 = ~all(centerImage == 0,2);
        squareT = double(square2)*double(square1);
        velR = regionprops(logical(squareT),'BoundingBox');
        centerImage = imcrop(centerImage,velR.BoundingBox);
        spaceChunk = imcrop(spaceChunk,velR.BoundingBox);
        timeChunk = imcrop(timeChunk,velR.BoundingBox);
        timeChunk_n = imcrop(timeChunk_n,velR.BoundingBox);

        

        sz1 = size(spaceChunk_vec);
        spaceChunkV = reshape(spaceChunk_vec,[prod(sz1(1:2)) sz1(3)]);

        sz1 = size(timeChunk_vec);
        timeChunkV = reshape(timeChunk_vec,[prod(sz1(1:2)) sz1(3)]);



        cur = [spaceChunkV timeChunkV];


        csTOT = PCA_REPROJ(cur,Etot,Utot);
        csTOT_view = reshape(csTOT,sz1);
       


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % color normalize
        for k = 1:size(csTOT_view,3)
            tmp = csTOT_view(:,:,k);
            tmp(maskIDX) = bindVec(tmp(maskIDX));
            %tmp(maskIDX) = tmp(maskIDX) - mean(tmp(maskIDX));
            csTOT_view(:,:,k) = tmp;
        end
        csTOT_view = imcrop(csTOT_view,velR.BoundingBox);

        if exist('gm')
            curZ = bsxfun(@minus,cur,masterU);
            curZ = bsxfun(@times,curZ,masterS.^-1);
            [kidx,nlogl,Pk] = gm.cluster(curZ);
            kidx = reshape(kidx,sz1(1:2));
            kidx = label2rgb(kidx);

            %Pk = reshape(Pk,sz1);
            %imshow(kidx,[]);
        end


        centerImage = repmat(centerImage,[1 1 3]);
        %imshow([centerImage spaceChunk timeChunk],[]);
        imshow([[centerImage spaceChunk];[csTOT_view timeChunk]],[]);
        drawnow



        masterSpace = [masterSpace ; cur];


    end
end

[masterSpaceZ,masterU,masterS] = zscore(masterSpace);
[Utot,Etot,Ltot] = PCA_FIT_FULLws(masterSpaceZ,3);

options = statset('Display','iter');
gm = fitgmdist(masterSpaceZ,5,'Options',options,'RegularizationValue',.01);
%%
%{
LSTM-network does not regress
%%
Y = totalCounts(mIDX);
for e = 1:size(v,3)
    X{e} = v(:,:,e);
end
%%
layers = [ ...
    sequenceInputLayer(size(v,1))
    lstmLayer(5,'OutputMode','last')
    fullyConnectedLayer(9)
    regressionLayer];

options = trainingOptions('sgdm', ...
    'MiniBatchSize',128,...
    'Shuffle','every-epoch',...
    'InitialLearnRate',.01,...
    'MaxEpochs',10, ...
    'ExecutionEnvironment','cpu',...
    'Plots','training-progress');

net = trainNetwork(X,Y,layers,options);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read the whole videos from disk
videoFile = '~/D1-3.hsv';
M = VideoReader(videoFile);
T = M.NumberofFrames;
H = M.Height;
W = M.Width;
S = zeros(H,W,T,'uint8');
% reset the video reader
cnt = 1;
M = VideoReader('~/D1-3.hsv');
while hasFrame(M)
    tmp = M.readFrame();
    S(:,:,cnt) = tmp(:,:,1);
    cnt = cnt + 1;
    cnt/T
end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
close all
for e = 1:size(S,3)
    imshow(S(:,:,e),[])
    title(num2str(e))
    drawnow
end
%% remove first N frames from stack
N = [45 1750];
[S] = removeFrontANDBack(S,N);
%% make plate mask
close all
plateMSK = zeros(size(S,1),size(S,2));
cp = hsize(plateMSK);

plateMSK(cp(1),cp(2)) = 1;
plateD = bwdist(plateMSK);
imshow(plateD,[]);
plateMSK = plateD < 500;
imshow(plateMSK,[]);
%% 
S = bsxfun(@times,double(S),plateMSK);
%% get the global stats
% each pixel falls within a range - normalize the pixel to be [0 1];
N = 10000;
uS = mean(single(S),3);
uSm = min(single(S),[],3);
uSM = max(single(S),[],3);
uRng = uSM - uSm;

sS = std(single(S),1,3);
dS = (diff(single(S),1,3));

udSm = min(single(dS),[],3);
udSM = max(single(dS),[],3);
udRng = udSM - udSm;


dS1 = zeros(size(S));
dS2 = zeros(size(S));
for e = 1:size(S,3)
    [dS1(:,:,e),dS2(:,:,e)] = gradient(single(S(:,:,e)));
    
    e
end

dTOT = (dS1.^2 + dS2.^2).^.5;

dTOTSm = min(single(dTOT),[],3);
dTOTSM = max(single(dTOT),[],3);
dTOTRng = dTOTSM - dTOTSm;

d1Sm = min(single(dS1),[],3);
d1SM = max(single(dS1),[],3);
d1Rng = d1SM - d1Sm;

d2Sm = min(single(dS2),[],3);
d2SM = max(single(dS2),[],3);
d2Rng = d2SM - d2Sm;

sV = std(dS(1:N));
uV = mean(dS(1:N));
%% distribution of window lengths
close all
vec = uRng(:);
vec = vec(randperm(numel(vec)));
vec = vec(1:100000);
options = statset('MaxIter',10);
window1Model = fitgmdist(vec,2,'Options',options);
[c1,~,pc1] = window1Model.cluster(uRng(:));
c1 = reshape(c1,size(uRng));
pc1 = reshape(pc1(:,2),size(uRng));
imshow(c1,[]);
close all
imshow(pc1,[]);
%% distribution of window widths
close all
vec = udRng(:);
vec = vec(randperm(numel(vec)));
vec = vec(1:10000);
ksdensity(vec);
window1Model = fitgmdist(vec,2,'Options',options);
c2 = window1Model.cluster(udRng(:));
c2 = reshape(c2,size(udRng));
imshow(c2,[]);
%% show both
close all
imshow((c1 == 2 & c2 == 2),[]);
%% sample places of interest and record x,y
X = [];
close all

b1 = linspace(0,1,50); % percent in intensity window
b2 = linspace(0,1,20); % percent in velocity window
b3 = linspace(0,1,30); % percet in gradient window
b4 = linspace(0,256,100); % abs location in intensity


statePDF = zeros(numel(b1),numel(b2),numel(b3),numel(b4));



X = [];
for e = 1:10:size(S,3)
    
    
    
N = 500000;
fidx = find((c1 == 2 & c2 == 2));
fidx = fidx(randperm(numel(fidx)));
fidx = fidx(1:N);

    
    q1 = single(S(:,:,e));
    q2 = single(dS(:,:,e));
    q3 = single(dTOT(:,:,e));
    q4 = uRng;
    
    
    s1 = (q1 - uSm).*uRng.^-1;
    s2 = (q2 - udSm).*udRng.^-1;
    s3 = (q3 - dTOTSm).*dTOTRng.^-1;
    s4 = q4;
    
    
    s1 = s1(fidx);
    s2 = s2(fidx);
    s3 = s3(fidx);
    s4 = s4(fidx);
    
    s1 = discretize(s1,b1);
    s2 = discretize(s2,b2);
    s3 = discretize(s3,b3);
    s4 = discretize(s4,b4);
    
    
    idx = sub2ind(size(statePDF),s1,s2,s3,s4);
    
    
    
    idx(isnan(idx)) = [];
    
    [sx,s2,s3,s4] = ind2sub(size(statePDF),idx);
    
    statePDF(idx) = statePDF(idx) + 1;
    
    %X = [X;[sx,sv,sg]];
    
    %vstatePDF = statePDF .* sum(statePDF(:))^-1;
    
    %vstatePDF = imfilter(vstatePDF,fspecial('disk',11),'replicate');
    
    %imshow(interp2(vstatePDF),[]);
    
    %plot(tmpX(:,1),tmpX(:,2),'.');
    %axis([0 1 0 1])
    drawnow
    
    %X = [X;tmpX];
    e
end
%%
B = ones(11,11,11,11);
B = B / sum(B(:));
%statePDFN = convn(statePDF,B,'same');
statePDFN = imfilter(statePDF,B,'replicate');
statePDFN = statePDFN * sum(statePDFN(:))^-1;
%%
statePDFN = statePDF * sum(statePDF(:))^-1;
%%
statePDFN = smooth3(statePDF,'gaussian',[21,21,21],11);
statePDFN = statePDFN * sum(statePDFN(:))^-1;


%{
%%
statePDFN = imfilter(statePDFN,fspecial('disk',11),'replicate');
statePDFN = statePDF * sum(statePDF(:))^-1;
%}
%{
%% find time trends
BX = bsxfun(@minus,S,uSm);
BX = bsxfun(@times,BX,uRng.^-1);

BV = bsxfun(@minus,dS,udSm);
BV = bsxfun(@times,BV,udRng.^-1);

sz1 = size(BX);
sz2 = size(BV);
BX = reshape(BX,[prod(sz1(1:2)) sz1(3)]);
BV = reshape(BV,[prod(sz2(1:2)) sz2(3)]);

kidx1 = ~any(isnan(BX),2);
kidx2 = ~any(isnan(BV),2);

uBX = mean(BX(kidx1,:),1);
uBV = mean(BV(kidx2,:),1);
%%
close all
plot(uBX);
figure;
plot(uBV);
%}
%% make background
warning off
filter_std_set = [5 5];
area_Filter = [20 2400];

perThreshDrop = .9;
diskSize = 11;

% prepare background(s)
BK = double(uS);
SS = imfilter(sS,fspecial('gaussian',[31 31],11),'replicate');

BK = imfilter(uSm,fspecial('gaussian',[31 31],3),'replicate');

%BK = imfilter(BK,fspecial('gaussian',[31 31],filter_std_set(1)),'replicate');
%BK = imopen(BK,strel('disk',diskSize,0));
%BK = imfilter(BK,fspecial('gaussian',[31 31],filter_std_set(2)),'replicate');


% play frame - BK
close all
velR = {};

disp = false;
disp2 = true;
clear R
h1 = figure;
h2 = figure;
p1 = [];
p2 = [];
interestMSK = (c1 == 2 & c2 == 2);
X = [];
for e = 1:(size(S,3)-1)
    
    q1 = single(S(:,:,e));
    q2 = single(dS(:,:,e));
    q3 = single(dTOT(:,:,e));
    q4 = uRng;
    
    
    s1 = (q1 - uSm).*uRng.^-1;
    s2 = (q2 - udSm).*udRng.^-1;
    s3 = (q3 - dTOTSm).*dTOTRng.^-1;
    s4 = q4;
    
    s1 = discretize(s1,b1);
    s2 = discretize(s2,b2);
    s3 = discretize(s3,b3);
    s4 = discretize(s4,b4);
    
    
    idx = sub2ind(size(statePDFN),s1,s2,s3,s4);
    
    
    
    idx(isnan(idx)) = 1;
    
    
    Prob = statePDFN(idx);
    Prob = reshape(Prob,size(q1));
    
    %Prob = Prob .* interestMSK;
    
    %Prob(find(interestMSK == 1)) = 1;
    
    mnV = max(Prob(find(interestMSK==1)));
    
    %mini = min(Prob(find(interestMSK==1 & Prob == 0)));
    %Prob(Prob == 0) = mini;
    
    Prob(interestMSK==0) = mnV;
    
    
    Prob = -log(Prob);
    
    
    %Prob(interestMSK==0) = mnV;
    
    
    %Prob = Prob.*interestMSK;
    
    %Prob(isnan(Prob(:))) = 0;
    %Prob(isinf(Prob(:))) = 1;
    
    
    %Prob = bindVec(Prob);
    
    
    
    
    fv = Prob(interestMSK==1);
    
    fv(isinf(fv)) = max(fv(~isinf(fv)));
    
    fv = bindVec(fv);
    th = graythresh(fv);
    th = .08;
    [py,pidx] = ksdensity(fv);
    [py,pidx] = ksdensity(fv,linspace(pidx(1),pidx(end),500));
    bin = imdilate(py,ones(1,21)) == py;
    bidx = find(bin);
    by = py(bidx);
    [~,bb] = sort(by);
    th = 1.5*pidx(bidx(bb(end-1)));
    
    th;
    
    fv = fv > th;
    msk = zeros(size(Prob));
    msk(interestMSK==1) = fv;
    msk = imerode(msk,strel('disk',1,0));
    
  
    
    [msk] = bwareaopenRange(logical(msk),[75 1500]);
    msk = imdilate(msk,strel('disk',1,0));
    velR = regionprops(logical(msk),'Centroid','Area','MajorAxisLength','PixelIdxList');
    %R = R([R.MajorAxisLength]<100);
    msk = zeros(size(msk));
    for o = 1:numel(velR)
        msk(velR(o).PixelIdxList) = 1;
    end
    
    %{
    out = flattenMaskOverlay(double(bindVec(q1)),logical(msk));
    %out = flattenMaskOverlay(out,logical(interestMSK),.8,'b');
    
    imshow(out,[]);
    hold on
    for o = 1:numel(R)
        plot(R(o).Centroid(1),R(o).Centroid(2),'m*');
        %text(R(o).Centroid(1),R(o).Centroid(2),num2str(R(o).MajorAxisLength),'Background','w');
    end
    
    hold off
    
    drawnow
    %}
    
    for o = 1:numel(velR)
        X = [X;[velR(o).Centroid e]];
    end
    e
    %waitforbuttonpress
end
%{

    
    % mean background subtraction
    tmp = double(S(:,:,e)) - BK;
    tmp = tmp .* (SS+1).^-1;
    
    
    
    
    vtmp = abs((dS(:,:,e) - uV)/sV);
    
    
    
    
    
    
    
    %tmp = imfilter(tmp,fspecial('gaussian',[31 31],11),'replicate');
    % zscore
   
    
    tmp = bindVec(tmp);
    msk = tmp > perThreshDrop*graythresh(tmp);
    vmsk = vtmp > 15;
    
    p1(e) = sum(vmsk(:))/numel(tmp);
    p2(e) = sum(msk(:) | vmsk(:))/numel(tmp);
    
    figure(h2);
    hold on
    plot(p1,'r')
    plot(p2,'b');
    hold off
    
    tmsk = msk | vmsk;
    
    
    
    midx1 = find(tmsk==1);
    midx0 = find(tmsk==0);
    midx0 = midx0(randperm(numel(midx0)));
    midx0 = midx0(1:numel(midx1));
    
    
    if disp
        subplot(1,2,1);
        plot(tmp(midx1),vtmp(midx1),'r.');
        hold on
        plot(tmp(midx0),vtmp(midx0),'b.');
        hold off 
        drawnow
    end

    
    R{e} = regionprops(logical(msk),'Centroid','Area','Perimeter');
    
    
    
    if disp2
        
        figure(h1);
        out = flattenMaskOverlay(bindVec(double(S(:,:,e))),msk,.5,'b');
        out = flattenMaskOverlay(out,vmsk,.5,'r');

        %subplot(1,2,2);

        imshow(out,[]);
        
        
        %imshow(bindVec(tmpv),[]);
        %drawnow
    end
    e
end
%}
%{
        %{
       
        
        
        %tmp = bindVec(tmp);
        msk = tmp > perThreshDrop*graythresh(tmp);
%}
        
        
        
        
        R = regionprops(msk,'Area','Perimeter','Centroid','PixelidxList');
        
        
        %fidx = find([R.Area] > area_Filter(1) & [R.Area] < area_Filter(2));
        %R = R(fidx);
        
        
        msk = zeros(size(msk));
        for o = 1:numel(fidx)
            msk(R(o).PixelIdxList) = 1;
        end

        
        for o = 1:numel(R)
            features{e}(:,o) = [R(o).Centroid';R(o).Area;R(o).Perimeter];
        end
        
        out = flattenMaskOverlay(S(:,:,e),logical(msk));

        imshow(out,[]);
        drawnow
    end
end
%}
%%
D = [X(:,3) X(:,1) X(:,2)];
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
velR = 20;
A = Radjacency(D(:,1:3)',velR);
%% 
A = sparse(size(D,1),size(D,1));
IDX = [];
V = [];
close all
bu = [];
for t = 1:(numel(UQ)-1)
    
    
    f1 = find(D(:,end) == UQ(t));
    f2 = find(D(:,end) == UQ(t+1));
    s1 = D(f1,1:3);
    s2 = D(f2,1:3);
    d = pdist2(s1,s2);
    
    b(t) = min(d(:));
    %d = tril(d);
    zidx = find(d~=0);
    g = d(zidx);
    [g,sidx] = sort(g);
    [idx1,idx2] = ind2sub(size(d),zidx(sidx(1:10)));
    bu(t) = mean(g(1:10));
    %{
    if t == 21
        figure;
        imshow(S(:,:,t),[]);
        hold on
        plot(s1(idx1,2),s1(idx1,3),'r.')
        plot(s2(idx2,2),s2(idx2,3),'go')
        plot(s2(:,2),s2(:,3),'g.')
        pause(1)
        figure;
        imshow(S(:,:,t+1),[]);
        hold on
        plot(s1(idx1,2),s1(idx1,3),'r.')
        plot(s2(idx2,2),s2(idx2,3),'go')
        waitforbuttonpress
        
    end
    %}
    
    d(d > 20) = 0;
    f = cartprod(f1,f2);
    %f = [f(:,2) f(:,1)];
    IDX = [IDX;f];
    V = [V;d(:)];
    t
end 
A = sparse(IDX(:,1),IDX(:,2),V,size(D,1),size(D,1));
%%
for e = 1:numel(UQ)
    fidx = find(D(:,end) == UQ(e));
    for i = 1:numel(fidx)
        for j = 1:numel(fidx)
             A(fidx(i),fidx(j)) = 0;
             A(fidx(j),fidx(i)) = 0;
        end
    end
    e
end
%%
velR = sparse(3,3);
velR(1,2) = 1;
velR(2,3) = 1;
[g1 , g2] = dijkstra(velR , 3 , 1);
       
        
%%
close all
fidxStr = find(D(:,end) == 1);
fidxStp = find(D(:,end) == max(UQ));
pathcost = [];
path = {};
close all
for str = 21:numel(fidxStr)
    parfor stp = 1:numel(fidxStp) 
        stp
        [path{str,stp} , pathcost(str,stp)] = dijkstra(A ,fidxStp(stp),fidxStr(str));
       
    
    
        if ~isempty(path{str,stp})
            for t = 1:10:size(S,3)
                D(path{str,stp}(t),1:4);
                imshow(S(:,:,t),[]);
                hold on
                plot(D(path{str,stp},2),D(path{str,stp},3),'m');
                plot(D(path{str,stp}(t),2),D(path{str,stp}(t),3),'g.');
                hold off
                drawnow
                TAU(t) = D(path{str,stp}(t),4);
                %waitforbuttonpress
                t
            end
            %}
            %close all
            plot(D(:,2),D(:,3),'b.','MarkerSize',1)
            hold on
            plot(D(fidxStr(str),2),D(fidxStr(str),3),'g*','MarkerSize',5)
            plot(D(fidxStp(stp),2),D(fidxStp(stp),3),'ro','MarkerSize',5)
            plot(D(path{str,stp} ,2),D(path{str,stp} ,3),'m','LineWidth',3)
            hold off
            drawnow
            pause(.4)
        else
            imshow(S(:,:,1),[]);hold on
            plot(D(fidxStr(str),2),D(fidxStr(str),3),'g*','MarkerSize',5);
            %figure;
            
            drawnow
            hold off
        end
        
        
        
    end
end
%% extract the X,Y,T
X = [];
for e = 1:numel(velR)
    for o = 1:numel(velR{e})
        X = [X ; [velR{e}(o).Centroid e o]];
    end
end

%% 
[~,BOX] = imcrop(S(:,:,1));
close all
%% make background
bSZ = 21;
filter_std_set = [5 5];
area_Filter = [20 2400];
perThreshDrop = .8;
diskSize = bSZ;

% prepare background(s)
BK = double(uS);
BK = imfilter(BK,fspecial('gaussian',[31 31],filter_std_set(1)),'replicate');
BK = imopen(BK,strel('disk',diskSize,0));
BK = imfilter(BK,fspecial('gaussian',[31 31],filter_std_set(2)),'replicate');

%%
clear J
jBK = double(imcrop(BK,BOX));
for e = 1:size(S,3)
    J(:,:,e) = double(imcrop(S(:,:,e),BOX)) - round(jBK);
    %J(:,:,e) = imfilter(J(:,:,e),fspecial('gaussian',[bSZ bSZ],3),'replicate');
    e
end
%%
clear S
%% try entropy filter
CH = 5;
FR = [];
N = 1200;
MSZ = numel(1:10:(size(J,3)-CH))*N;
%H = zeros(2*256,2*256,MSZ);
yIDX = 200:300;
xIDX = 200:300;
H = zeros(numel(yIDX),numel(xIDX),MSZ);
cnt = 1;
bSZ = 51;
selectedTrainIDX = 1:(size(J,1)*size(J,2));
selectedTrainIDX = reshape(selectedTrainIDX,[size(J,1) size(J,2)]);
selectedTrainIDX = im2colF(selectedTrainIDX,[bSZ bSZ],[1 1]);
flag = true;



for e = 1:10:(size(J,3)-CH)
    e
    slice = J(:,:,e:e+CH);
    dSlice = diff(slice,1,3);
    STACK = [];

    %tmp = im2colF(double(slice(:,:,1)),[bSZ bSZ],[1 1]);
    ridx = randperm(size(selectedTrainIDX,2));

    tH = zeros(2*256,2*256,N);

    for s1 = 1:size(dSlice,3)
        tmp = double(slice(:,:,s1));
        dtmp = double(dSlice(:,:,s1));

        %tmp = im2colF(double(slice(:,:,s)),[bSZ bSZ],[1 1]);
        %dtmp = im2colF(double(dSlice(:,:,s)),[bSZ bSZ],[1 1]);
        %ridx = randperm(size(tmp,2));
        %tmp = tmp(:,ridx);
        %dtmp = dtmp(:,ridx);

        for v = 1:N

            tIDX = selectedTrainIDX(:,ridx(v));
            %ttH = zeros(256,2*256);
            %HH = tH(:,:,v);
            %ind = sub2ind([size(ttH,1) size(ttH,2)],tmp(:,v)+1,dtmp(:,v)+256);
            ind = sub2ind([size(tH,1) size(tH,2) size(tH,3)],tmp(tIDX)+256,dtmp(tIDX)+256,v*ones(numel(tIDX),1));
            for h = 1:numel(ind)
                %ttH(ind(h)) = ttH(ind(h)) + 1;
                %HH(ind(h)) = HH(ind(h)) + 1;
                tH(ind(h)) = tH(ind(h)) + 1;
            end
            %{
            HH = tH(:,:,v);
            ind = sub2ind([size(H,1) size(H,2)],tmp(:,v)+1,dtmp(:,v)+256);
            for h = 1:numel(ind)
                HH(ind(h)) = HH(ind(h)) + 1;
            end
            tH(:,:,v) = HH;
            %}
            %imshow(HH,[]);
            %drawnow
        end
        %STACK = [STACK;tmp];
    end

    
    
    %{
    for s = 1:N
        imshow(tH(:,:,s),[]);
        drawnow
    end
    %}
    
    H(:,:,cnt:(cnt+N-1)) = tH(yIDX,xIDX,:);
    cnt = cnt + N;
    
    %FR(:,e) = std(STACK,1,1);
end
%%
H(:,:,cnt:end) = [];

%% trim H
xIDX = 200:300;
a = squeeze(mean(sum(H,1),3));
b = squeeze(mean(sum(H,2),3));
plot(a);
%%
close all
yIDX = 200:300;
b = squeeze(mean(sum(H,2),3));
plot(b)
%%
%subH = H(yIDX,xIDX,:);
subH = H;
sz = size(subH);
subH = reshape(subH,[prod(sz(1:2)) sz(3)]);
subH = subH * bSZ.^-2;
%%
[U,E,L] = PCA_FIT_FULL_Tws(subH,3);
%%
C = PCA_REPROJ_T(subH,E,U);
%% sweep
close all
SW = sweepPCA(C',E,U',L.^.5,1,5);
SW = squeeze(SW);
for e = 1:size(SW,1)
    tmp = reshape(SW(e,:),[sz(1:2)]);
    %tmp = reshape(U,[sz(1:2)]);
    %tmp = reshape(SW(e,:)-U',[sz(1:2)]);
    imshow(tmp,[]);
    %waitforbuttonpress
    plot(mean(tmp,2))
    %waitforbuttonpress
end
%%
close all
rE = reshape(bsxfun(@plus,E,U),[sz(1:2) 3]);
imshow(rE(:,:,:),[]);
%%
NG = 3;
options = statset('Display','iter','MaxIter',5);
gmm = fitgmdist(C',NG,'Options',options,'RegularizationValue',.01,'Replicates',10);
%% try entropy filter
CH = 5;
MAS = [];
MAS2 = [];
for e = 1:(size(J,3)-CH)
    slice = J(:,:,e:e+CH);
    dSlice = diff(slice,1,3);
    STACK = [];
    
    tmp = im2colF(double(slice(:,:,1)),[bSZ bSZ],[1 1]);
    tH = zeros(2*256,2*256,size(tmp,2));
    
    for s1 = 1:size(dSlice,3)
        
        %tmp = double(slice(:,:,s));
        %dtmp = double(dSlice(:,:,s));

        
        tmp = im2colF(double(slice(:,:,s1)),[bSZ bSZ],[1 1]);
        dtmp = im2colF(double(dSlice(:,:,s1)),[bSZ bSZ],[1 1]);
        tic
        for v = 1:size(selectedTrainIDX,2)
            %tIDX = mIDX(:,v);
            %tic
            %ttH = zeros(256,2*256);
            %HH = tH(:,:,v);
            %ind = sub2ind([size(tH,1) size(tH,2)],tmp(:,v)+256,dtmp(:,v)+256);
            %ind = sub2ind([size(tH,1) size(tH,2) size(tH,3)],tmp(mIDX(:,v))+256,dtmp(mIDX(:,v))+256,v*ones(size(tIDX,1),1));
            ind = sub2ind([size(tH,1) size(tH,2) size(tH,3)],tmp(:,v)+256,dtmp(:,v)+256,v*ones(size(tmp,1),1));
            
            %uidx = unique(ind);
            %su = sum(ind==uidx');
            %tH(uidx) = tH(uidx) + su';
            
            for h = 1:numel(ind)
                %ttH(ind(h)) = ttH(ind(h)) + 1;
                %HH(ind(h)) = HH(ind(h)) + 1;
                tH(ind(h)) = tH(ind(h)) + 1;
            end
            
            %tH(:,:,v) = tH(:,:,v) + ttH;
            %tH(:,:,v) = tH(:,:,v) + ttH;
            %tH(:,:,v) = HH;
            %v
            %imshow(HH,[]);
            %drawnow
            %toc
        end
        toc
        s1
        %STACK = [STACK;tmp];
    end
    
    sub_tH = tH(yIDX,xIDX,:);
    sz = size(sub_tH);
    sub_tH = reshape(sub_tH,[prod(sz(1:2)) sz(3)]);
    tC = PCA_REPROJ_T(sub_tH,E,U);
    
    [kidx,nlogl,P] = gmm.cluster(tC');
    
    kidx = col2im(kidx,[bSZ bSZ],[size(slice,1) size(slice,2)]);
    KP = [];
    for k = 1:NG
        KP(:,:,k) = col2im(P(:,k),[bSZ bSZ],[size(slice,1) size(slice,2)]);
        %KP(:,:,k) = bindVec(K(:,:,k));
    end
    
    for k = 1:3
        K(:,:,k) = col2im(tC(k,:),[bSZ bSZ],[size(slice,1) size(slice,2)]);
        %K(:,:,k) = bindVec(K(:,:,k));
    end
    
    MAS(:,:,:,e) = KP;
    MAS2(:,:,:,e) = K;
    e
    cnt
    %FR(:,e) = std(STACK,1,1);
end
%%
close all
h1 = figure;
h2 = figure;
for e = 1:size(MAS,4)
   
    [~,midx] = max(MAS(:,:,:,e),[],3);
    msk = midx == 3;
   
    figure(h2)
    tmp = J(:,:,e);
    
    for r = 1:4
        tmp(1:25,:) = [];
        tmp = imrotate(tmp,90);
    end
    imshow(tmp,[])
    
    out = flattenMaskOverlay(bindVec(tmp),msk);
    figure(h1)
    imshow(out,[]);
    drawnow
    %waitforbuttonpress
end
%% get the mean
uS = mean(single(S),3);
%% match features
for e = 1:numel(features)
    f1 = features{e};
    f2 = features{e+1};
    D = zeros(size(f1,2),size(f2,2));
    for i = 1:size(f1,2)
        for j = 1:size(f2,2)
            delta = f1(:,i) - f2(:,j);
            D(i,j) = sum(delta(1:2).*delta(1:2),1).^.5;
        end
    end
    
    for i = 1:size(f1,2)
        [r,sidx] = sort(D(i,:));
        fgamma{e}(i,:) = sidx(1:3);
    end
    e
end
%% plot gamma
P = [];
for e = 1:size(S,3)
    if e == 1
        for pt = 1:size(features{e},2)
            P(pt,:,e) = features{e}(1:2,pt)';
        end
        
    else
        
        for pt = 1:size(P,1)
            P(pt,:,e) = features{e}(1:2,fgamma{e-1}(pt,1))';
        end
    end
    
    
    imshow(S(:,:,e),[]);
    hold on
    for pt = 1:size(P,1)
        trace = squeeze(P(pt,:,:));
        %plot(trace(1,:),trace(2,:),'r')
    end
    
    e
    
end
%% remove the mean frame from the stack
sS = std(single(S),1,3);
Zthreshold = 2;
for e = 1:size(S,3)
    MSK = (single(S(:,:,e)) - uS) > Zthreshold*sS;
    %S1(:,:,e) = double(MSK).*double(S(:,:,e)) + double(~MSK).*uS;
    S1(:,:,e) = double(MSK).*double(S(:,:,e));
    %imshow(S1(:,:,e),[]);
    %title(num2str(e))
    %drawnow
    %imshow(S1(:,:,e),[]);
    %drawnow
    e
    size(S,3)
end
%% make dt and ds
%S1 = bsxfun(@minus,single(S),uS);
DT = diff(S1,1,3);
for e = 1:size(S1,3)
    [d1,d2] = gradient(S1(:,:,e));
    DS(:,:,e) = (d1.^2 + d2.^2).^.5;
    e
end
%%
%{
%%
for e = 1:10%size(S,3)
    tmp = imfilter(S(:,:,e),fspecial('gaussian',[31 31],7),'replicate');
    sz = size(tmp);
    tmp = imresize(tmp,.25);
    tmp = imopen(tmp,strel('disk',7,0));
    tmp = imresize(tmp,sz);
    tmp = imfilter(tmp,fspecial('gaussian',[31 31],7),'replicate');
    tmp = S(:,:,e) - tmp;
    
    pS(:,:,e) = tmp;
    imshow(pS(:,:,e),[]);
    drawnow
end
%%
tmp = double(mean(S,3));
tmp = imfilter(tmp,fspecial('gaussian',[31 31],7),'replicate');
sz = size(tmp);
tmp = imresize(tmp,.25);
tmp = imopen(tmp,strel('disk',7,0));
tmp = imresize(tmp,sz);
tmp = imfilter(tmp,fspecial('gaussian',[31 31],7),'replicate');
BK = tmp;
for e = 1:10%size(S,3)
    pS(:,:,e) = double(S(:,:,e)) - BK;
    imshow(pS(:,:,e),[]);
    drawnow
end
%%
sD = diff(S,1,3);
%%
S = S(:,:,1:end-1);
%}
%% gather the feature space over some of the image patches
% image patches are drawn from S1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
patchSKUIP = 4;
W = [];
for e = 1:patchSKUIP:size(S1,3)
    tic
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for raw data
    tmp = S1(:,:,e);
    H0 = im2colF(double(tmp),[51 51],[11 11]);
    idx = randperm(size(H,2));
    H0 = H0(:,idx);
    H0 = H0(:,1:5000);
    H0 = sort(H0,1);
    H0 = H0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[tmp] = imcrop(sD(:,:,e),BOX);
    % for diff data
    tmp = DT(:,:,e);
    H1 = im2colF(double(tmp),[51 51],[11 11]);
    %idx = randperm(size(H,2));
    H1 = H1(:,idx);
    H1 = H1(:,1:5000);
    H1 = sort(H1,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[tmp] = imcrop(sD(:,:,e),BOX);
    % for diff data
    tmp = DS(:,:,e);
    H2 = im2colF(double(tmp),[51 51],[11 11]);
    %idx = randperm(size(H,2));
    H2 = H2(:,idx);
    H2 = H2(:,1:5000);
    H2 = sort(H2,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % cat the feature vectors together
    W = [W [H0;H1;H2]];
    
    %W = [W H];
    toc
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PCA on the data each seperate feature space - [X,dS,dT];
nD = 3;
%nD = 1;
CHUNK = size(W,1)/nD;
IDX1 = 1:CHUNK;
IDX2 = (CHUNK+1):(2*CHUNK);
IDX3 = (2*CHUNK+1):(3*CHUNK);
[U1,E1,L1] = PCA_FIT_FULL_Tws(W(IDX1,:),3);
[U2,E2,L2] = PCA_FIT_FULL_Tws(W(IDX2,:),3);
[U3,E3,L3] = PCA_FIT_FULL_Tws(W(IDX3,:),3);
L1 = L1(1).^.5;
L2 = L2(1).^.5;
L3 = L3(1).^.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% back project the data
Sc1 = PCA_REPROJ_T(W(IDX1,:),E1,U1);
Sc2 = PCA_REPROJ_T(W(IDX2,:),E2,U2);
Sc3 = PCA_REPROJ_T(W(IDX3,:),E3,U3);
Sc1 = bsxfun(@times,Sc1,L1^-1);
Sc2 = bsxfun(@times,Sc2,L2^-1);
Sc3 = bsxfun(@times,Sc3,L3^-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% direct sum the data together and perform PCA across spaces
ScT = [Sc1;Sc2;Sc3];
[UT,ET,LT] = PCA_FIT_FULL_Tws(ScT,size(ScT,1));
SccT = PCA_REPROJ_T(ScT,ET,UT);
close all
plot(cumsum(LT)/sum(LT));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% look at the PCA of the raw data - space together
sel = nchoosek(1:6,2);
for e = 1:size(sel,1)
    plot(SccT(sel(e,1),:),SccT(sel(e,2),:),'.');
    drawnow
    waitforbuttonpress
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% look at the PCA of the raw data - each space seperate
close all
figure
plot(Sc1(1,:),Sc1(2,:),'.');
figure
plot(Sc1(1,:),Sc1(3,:),'.');
figure
plot(Sc2(1,:),Sc2(2,:),'.');
figure
plot(Sc3(1,:),Sc3(2,:),'.');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% "straighten" the data
funcS{1} = @(value)min(value);
funcS{2} = @(value)max(value);
funcS{3} = @(value)mean(value);
funcS{4} = @(value)std(value);
funcS{5} = @(value).5*(max(value) - min(value)) + min(value);

% create combine function of the functions run on the data
combineFunction = @(data,resultY)(data - resultY(5,:));
combineFunction = @(data,resultY)(data - resultY(1,:));

%%%%%%%%%%%%%%%%
% make PC2 a function of PC1
% Data = X space
NSP = 200;
[funcY1] = straightMiron(Sc1,[1 2],funcS,NSP);
nSc1 = applyStraightenMiron(Sc1,funcY1,combineFunction);
close all
plot(Sc1(1,:),Sc1(2,:),'.')
hold on
plot(nSc1(1,:),nSc1(2,:),'r.')

%%%%%%%%%%%%%%%%
% make PC2 a function of PC1
% Data = dS space
NSP = 200;
[funcY2] = straightMiron(Sc2,[1 2],funcS,NSP);
nSc2 = applyStraightenMiron(Sc2,funcY2,combineFunction);
figure
plot(Sc2(1,:),Sc2(2,:),'.')
hold on
plot(nSc2(1,:),nSc2(2,:),'r.')


%%%%%%%%%%%%%%%%
% make PC2 a function of PC1
% Data = dT space
NSP = 200;
[funcY3] = straightMiron(Sc3,[1 2],funcS,NSP);
nSc3 = applyStraightenMiron(Sc3,funcY3,combineFunction);
figure
plot(Sc3(1,:),Sc3(2,:),'.')
hold on
plot(nSc3(1,:),nSc3(2,:),'r.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% investigate the relationship between PC1 and PC2
K = convhulln(Sc1');
%{
sweepNP = 100;
sweepDIM = 1;
sweepDomain = linspace(min(dataToSweep(sweepDIM,:)),max(dataToSweep(sweepDIM,:)),sweepNP);
sweepCurve(1:2,:) = [sweepDomain;funcY1(1).func(sweepDomain)];
%}
%{
close all
plot(dataToSweep(1,:),dataToSweep(2,:),'.')
hold on
plot(sweepCurve(1,:),sweepCurve(2,:),'r')
waitforbuttonpress
sweepF = PCA_BKPROJ_T(sweepCurve,E1,U1);
close all
plot(sweepF)
%}

close all
figure
plot(Sc1(1,:),Sc1(3,:),'.')
figure
plot(Sc1(1,:),Sc1(2,:),'.')
hold on
plot(Sc1(1,:),nSc1(2,:),'.')
figure
plot(Sc1(2,:),Sc1(3,:),'.')
hold on
plot(nSc1(2,:),Sc1(3,:),'.')
figure;
plot3(Sc1(1,:),Sc1(2,:),Sc1(3,:),'.')
hold on;
trisurf(K,Sc1(1,:)',Sc1(2,:)',Sc1(3,:)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sweep a select line in PC space to see what it does to the eigen function
%% first sweep the orginal - PC1
dataToSweep = Sc1;
sweepDIM = 1;
sweepNP = 10;
sweepRange = std(Sc1,1,2);
sweepDomain = linspace(-sweepRange(sweepDIM),sweepRange(sweepDIM),sweepNP);
sweepDomain = linspace(min(dataToSweep(sweepDIM,:)),max(dataToSweep(sweepDIM,:)),sweepNP);
sweepCurve = zeros(size(dataToSweep,1),sweepNP);
sweepCurve(sweepDIM,:) = L1*sweepDomain;
sweepF = PCA_BKPROJ_T(sweepCurve,E1,U1);
plot(sweepF)
waitforbuttonpress
%% first sweep the orginal - PC2
dataToSweep = Sc1;
sweepDIM = 2;
sweepNP = 10;
sweepRange = std(Sc1,1,2);
sweepDomain = linspace(-sweepRange(sweepDIM),sweepRange(sweepDIM),sweepNP);
sweepDomain = linspace(min(dataToSweep(sweepDIM,:)),max(dataToSweep(sweepDIM,:)),sweepNP);
sweepCurve = zeros(size(dataToSweep,1),sweepNP);
sweepCurve(sweepDIM,:) = L1*sweepDomain;
sweepF = PCA_BKPROJ_T(sweepCurve,E1,U1);
plot(sweepF)
waitforbuttonpress
%% look at the min curve
sweepDIM = 1;
sweepDomain = linspace(min(dataToSweep(sweepDIM,:)),max(dataToSweep(sweepDIM,:)),sweepNP);
sweepCurve(1:2,:) = [sweepDomain;funcY1(1).func(sweepDomain)];
close all
plot(dataToSweep(1,:),dataToSweep(2,:),'.')
hold on
plot(sweepCurve(1,:),sweepCurve(2,:),'r')
waitforbuttonpress
sweepF = PCA_BKPROJ_T(sweepCurve,E1,U1);
close all
plot(sweepF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% look at the crosses between [X,dS,dT]
close all
plot(nSc1(1,:),nSc3(1,:),'.')
figure
plot(nSc3(1,:),nSc1(1,:),'.')
figure;
plot(nSc1(2,:),nSc3(2,:),'.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% twist the crosses - make the PC's of an element of [X,dS,dT] a function of the other
close all
NSP = 500;
dataX1 = cat(1,nSc3(1,:),nSc1(1,:));
[funcY3_cross1] = straightMiron(dataX1,[1 2],funcS,NSP);
ndataX1 = applyStraightenMiron(dataX1,funcY3_cross1,combineFunction);
[funcY3_cross2] = straightMiron(ndataX1,[2 1],funcS,NSP);
ndataX2 = applyStraightenMiron(ndataX1,funcY3_cross2,combineFunction);
figure
plot(dataX1(1,:),dataX1(2,:),'.')
hold on
figure
plot(ndataX1(1,:),ndataX1(2,:),'r.')
figure;
plot(ndataX2(1,:),ndataX2(2,:),'g.')
[Z] = myPDF2(ndataX1,[100 100]);
figure;
lZ = -log(Z);
lZ(isinf(lZ(:))) = 0;
imshow(lZ,[]);
%{
%%%%%%%%%%%%%%%%%%%%%%
% NOTE THAT GMM @ THIS POINT SEEM TO BE A BIT EARLY - 
% RATHER NEEDING TO EXPLORE - AS OF SEP 25,  2019
%%%%%%%%%%%%%%%%%%%%%%
%% fit GMM model(s) to the data
close all
NG = 2;
dataToFit = ndataX1;
kidx = kmeans(dataToFit',NG);

options = statset('Display','iter','MaxIter',1000);
gm = fitgmdist(dataToFit',2,'Options',options,'Replicates',2,'Start','plus');
gmPDF = @(x,y)pdf(gm,[x y]);
figure;
CL = {'r.','g.','b.'};
hold on
for k = 1:NG
    %scatter(dataToFit(1,kidx==k),dataToFit(2,kidx==k),10,CL{k})
end
kidx2 = gm.cluster(dataToFit');
hold on
for k = 1:NG
    scatter(dataToFit(1,kidx2==k),dataToFit(2,kidx2==k),10,CL{k})
end
[n1,n2] = ndgrid(linspace(min(dataToFit(1,:)),max(dataToFit(1,:)),100),linspace(min(dataToFit(2,:)),max(dataToFit(2,:)),100));
E = gmPDF(n1(:),n2(:));
E = reshape(E,size(n1));
hold on
contour(gca,n1,n2,E,10);
hold on
idxL = kidx2;
BKIDX = mode(idxL);
%}
%% second level PCA
subIDX = idxL~=BKIDX;
CHUNK = size(W,1)/3;
IDX1 = 1:CHUNK;
IDX2 = (CHUNK+1):(2*CHUNK);
IDX3 = (2*CHUNK+1):(3*CHUNK);
[U1_2,E1_2,L1_2] = PCA_FIT_FULL_Tws(W(IDX1,subIDX),3);
[U2_2,E2_2,L2_2] = PCA_FIT_FULL_Tws(W(IDX2,subIDX),3);
[U3_2,E3_2,L3_2] = PCA_FIT_FULL_Tws(W(IDX3,subIDX),3);
L1_2 = L1_2(1);
L2_2 = L2_2(1);
L3_2 = L3_2(1);
%% back project the data
Sc1_2 = PCA_REPROJ_T(W(IDX1,subIDX),E1_2,U1_2);
Sc2_2 = PCA_REPROJ_T(W(IDX2,subIDX),E2_2,U2_2);
Sc3_2 = PCA_REPROJ_T(W(IDX3,subIDX),E3_2,U3_2);
Sc1_2 = bsxfun(@times,Sc1_2,L1_2^-1);
Sc2_2 = bsxfun(@times,Sc2_2,L2_2^-1);
Sc3_2 = bsxfun(@times,Sc3_2,L3_2^-1);
%% look at the PCA of the raw data
close all
figure
plot(Sc1_2(1,:),Sc1_2(2,:),'.');
figure
plot(Sc2_2(1,:),Sc2_2(2,:),'.');
figure
plot(Sc3_2(1,:),Sc3_2(2,:),'.');
%% "straighten" the data after section out the data
funcS{1} = @(value)min(value);
funcS{2} = @(value)max(value);
funcS{3} = @(value)mean(value);
funcS{4} = @(value)std(value);
combineFunction = @(data,resultY)(data - resultY(1,:));

NSP = 200;
[funcY1_2] = straightMiron(Sc1_2,[1 2],funcS,NSP);
nSc1_2 = applyStraightenMiron(Sc1_2,funcY1_2,combineFunction);
close all
plot(Sc1_2(1,:),Sc1_2(2,:),'.')
hold on
plot(nSc1_2(1,:),nSc1_2(2,:),'r.')


NSP = 200;
[funcY2_2] = straightMiron(Sc2_2,[1 2],funcS,NSP);
nSc2_2 = applyStraightenMiron(Sc2_2,funcY2_2,combineFunction);
figure
plot(Sc2_2(1,:),Sc2_2(2,:),'.')
hold on
plot(nSc2_2(1,:),nSc2_2(2,:),'r.')



NSP = 200;
[funcY3_2] = straightMiron(Sc3_2,[1 2],funcS,NSP);
nSc3_2 = applyStraightenMiron(Sc3_2,funcY3_2,combineFunction);
figure
plot(Sc3_2(1,:),Sc3_2(2,:),'.')
hold on
plot(nSc3_2(1,:),nSc3_2(2,:),'r.')

%% look at the crosses
close all
plot(nSc1_2(1,:),nSc3_2(1,:),'.')
figure
plot(nSc3_2(1,:),nSc1_2(1,:),'.')
figure;
plot(nSc1_2(2,:),nSc3_2(2,:),'.')
%% twist the crosses
close all
NSP = 500;
dataX1 = cat(1,nSc3_2(1,:),nSc1_2(1,:));
[funcY3_cross1] = straightMiron(dataX1,[1 2],funcS,NSP);
ndataX1 = applyStraightenMiron(dataX1,funcY3_cross1,combineFunction);
[funcY3_cross2] = straightMiron(ndataX1,[2 1],funcS,NSP);
ndataX2 = applyStraightenMiron(ndataX1,funcY3_cross2,combineFunction);
figure
plot(dataX1(1,:),dataX1(2,:),'.')
hold on
figure
plot(ndataX1(1,:),ndataX1(2,:),'r.')
figure;
plot(ndataX2(1,:),ndataX2(2,:),'g.')
[Z] = myPDF2(ndataX1,[100 100]);
figure;
lZ = -log(Z);
lZ(isinf(lZ(:))) = 0;
imshow(lZ,[]);

%% fit GMM model(s) to the data
close all
NG = 2;
dataToFit = ndataX1;
kidx = kmeans(dataToFit',NG);

options = statset('Display','iter','MaxIter',1000);
gm = fitgmdist(dataToFit',2,'Options',options,'Replicates',2,'Start','plus');
gmPDF = @(x,y)pdf(gm,[x y]);
figure;
CL = {'r.','g.','b.'};
hold on
for k = 1:NG
    %scatter(dataToFit(1,kidx==k),dataToFit(2,kidx==k),10,CL{k})
end
kidx2 = gm.cluster(dataToFit');
hold on
for k = 1:NG
    scatter(dataToFit(1,kidx2==k),dataToFit(2,kidx2==k),10,CL{k})
end
[n1,n2] = ndgrid(linspace(min(dataToFit(1,:)),max(dataToFit(1,:)),100),linspace(min(dataToFit(2,:)),max(dataToFit(2,:)),100));
E = gmPDF(n1(:),n2(:));
E = reshape(E,size(n1));
hold on
contour(gca,n1,n2,E,10);
hold on
idxL = kidx2;
BKIDX = mode(idxL);

%% look at scatter plot of sc1
close all
plot(Sc1_2(1,:),Sc1_2(2,:),'.')
figure;
close all
plot(Sc1_2(2,:),Sc2_2(2,:),'.')
%%
%% second level cluster
close all
grps2 = 4;
GMModel_2 = fitgmdist([Sc1_2;Sc2_2]',grps2);
[idxL_2] = cluster(GMModel_2,[Sc1_2;Sc2_2]');
CL = {'r.','g.','b.','k.'};
h1 = figure;
h2 = figure;
for uTrain = 1:grps2
    sidx = idxL_2 == uTrain;
    figure(h1);
    plot(Sc1_2(1,sidx),Sc2_2(1,sidx),CL{uTrain});
    hold on
    figure(h2)
    ksdensity(Sc1_2(1,sidx))
    hold on
end
%%
UQ = unique(idxL);
selIDX = idxL ~= BKIDX;
GMModel2 = fitgmdist([Sc1(:,selIDX);Sc2(:,selIDX)]',3);
%%
plot(sum(W(IDX1,:),1),Sc1(1,:),'.');
plot(Sc1(1,:),Sc1(2,:),'.');
%%
ds = 5;
cnt = 1;
vec1 = [];
STORE = [];
TARN = [];
NN = [];
cnt = 1;
numIMG = size(S1,3)
numIMG = 100;
%numIMG = 10;
velR = 15;
N = 25;
MSZ = 500;
MASK = zeros(size(S1,1),size(S1,2));
MASK(end/2,end/2) = 1;
MASK = bwdist(MASK) < MSZ;
for e = 1:1:numIMG
    tic
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for raw data
    tmp = S1(:,:,e);
    reOLD = size(tmp);
    
    NN(e) = randi([15,30],1,1);
    tmp = generateImageP(MASK,velR,NN(e));
    imshow(tmp,[]);
    drawnow
    H1 = im2colF(double(tmp),[51 51],[ds ds]);
    
    
    
    
    %H1 = im2colF(double(tmp),[51 51],[(51-1) (51-1)]);
    H1 = sort(H1,1);
    
    % generte targets for detectors
    TARN(cnt,:) = sum(H1,1)*(pi*velR)^-2;
    
    
    
    Sc1 = PCA_REPROJ_T(H1,E1,U1);
    Sc1 = bsxfun(@times,Sc1,L1.^-1);
    
    %Sc1_2 = PCA_REPROJ_T(H1,E1_2,U1_2);
    %Sc1_2 = bsxfun(@times,Sc1_2,L1_2.^-1);
    
    STORE(:,:,cnt) = Sc1;
    cnt = cnt + 1;
    
    %[Count,img] = countOP(Sc1,x,25);
    % out = 
    
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[tmp] = imcrop(sD(:,:,e),BOX);
    % for diff data
    tmp = D1(:,:,e);
    H2 = im2colF(double(tmp),[51 51],[ds ds]);
    H2 = sort(H2,1);
    
    Sc2 = PCA_REPROJ_T(H2,E2,U2);
    Sc2 = bsxfun(@times,Sc2,L2.^-1);
    
    Sc2_2 = PCA_REPROJ_T(H2,E2_2,U2_2);
    Sc2_2 = bsxfun(@times,Sc2_2,L2_2.^-1);
   
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % diff data
    tmp = DS(:,:,e);
    reOLD = size(tmp);
    H3 = im2colF(double(tmp),[51 51],[ds ds]);
    H3 = sort(H3,1);
    
    Sc3 = PCA_REPROJ_T(H3,E3,U3);
    Sc3 = bsxfun(@times,Sc3,L3.^-1);
    Sc3_2 = PCA_REPROJ_T(H3,E3_2,U3_2);
    Sc3_2 = bsxfun(@times,Sc3_2,L3_2.^-1);
    
    vec2(cnt,:) = mean(Sc3,2)';
    cnt = cnt + 1;
    
    
    
    %[idx1,~,P1] = cluster(GMModel,[Sc1;Sc2;Sc3]');
    [idx1,~,P1] = cluster(GMModel,[Sc1]');
    [idx2,~,P2] = cluster(GMModel_2,[Sc1_2;Sc2_2]');
   
    for k = 1:size(P1,2)
        IMG1(:,:,k,e) = col2im(P1(:,k),[51 51],size(tmp));
        IMG_LAB1(:,:,e) = col2im(idx1,[51 51],size(tmp));
    end
    
    
    for k = 1:size(P2,2)
        IMG2(:,:,k,e) = col2im(P2(:,k),[51 51],size(tmp));
        IMG_LAB2(:,:,e) = col2im(idx2,[51 51],size(tmp));
    end
    
    vec1(cnt,:) = mean(Sc1,2)';
    %}
    
    %{
    % back apply beta
    PM = beta'*[Sc1;Sc3];
    PM = col2im(PM,[51 51],size(tmp));
    %}
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[tmp] = imcrop(sD(:,:,e),BOX);
    % for diff data
    tmp = D1(:,:,e);
    H = im2colF(double(tmp),[51 51],[ds ds]);
    H = sort(H,1);
    Sc2 = PCA_REPROJ_T(H,E2,U2);
    Sc2 = bsxfun(@times,Sc2,L2.^-.5);
   
    
    WSc = [Sc1;Sc2];
    [idx,nlogl,P] = cluster(GMModel,WSc');
    
    [idx2,nlogl2,P2] = cluster(GMModel2,WSc');
    
    
    
    LAB1 = col2im(idx,[51 51],size(tmp));
    for k = 1:size(P,2)
        tmp1(:,:,k) = col2im(P(:,k),[51 51],reOLD);
    end
    
    
    LAB2 = col2im(idx2,[51 51],size(tmp));
    for k = 1:size(P2,2)
        tmp2(:,:,k) = col2im(P2(:,k),[51 51],reOLD);
    end
    
    
    
    IMG(:,:,:,e) = bsxfun(@times,tmp2,LAB1~=BKIDX);
    
    
    
    imshow(IMG(:,:,:,e),[]);
    drawnow
    toc
    %}
    toc
end
%%
func = @(vec)countOP(STORE,vec,25);
options = optimoptions('fmincon','Display','iter');
x0 = 2*(rand(1,5)-.5);
x0(end) = 10000/2;
xTrain = fmincon(func,x0,[],[],[],[],[-1;-1;-1;-100;10^3],[1;1;1;100;10^4],[],options);
%%
options = optimset('Display','iter');
xTrain = fminsearch(func,2*(ones(1,5)-.5),options);
%%
func = @(vec)countOP(STORE,vec,25);
options = optimoptions('particleswarm','Display','iter');
xTrain = particleswarm(func,5,[-1;-1;-1;-100;10^3],[1;1;1;100;10^4],options);
%%
func = @(vec)countOP(STORE,vec,25);
options = optimoptions('particleswarm','Display','iter');
xTrain = particleswarm(func,5,[-1;-1;-1;-100;10^3],[1;1;1;100;10^4],options);
%%
sz = [3 3];
func = @(vec)countOP3(STORE,vec,25,sz);
options = optimoptions('particleswarm','Display','iter','MaxIter',200);
xo = [rand(1,prod(sz)) 0 0 0  10 10 10];
LB = [-ones(prod(sz),1);-10*ones(sz(1),1);.01*ones(sz(1),1)];
UB = [ones(prod(sz),1);10*ones(sz(1),1);10*ones(sz(1),1)];
xTrain = particleswarm(func,numel(xo),LB,UB,options);
%%
global BD
BD = [];
sz = [3 3];
func = @(vec)countOP4(STORE,vec,NN,TARN,[],sz);
OF = @(A,B)outputPLT(func,A,B);
options = optimoptions('particleswarm','Display','iter','MaxIter',200,'UseParallel',true,'PlotFcn',OF);
UU = mean(mean(STORE,2),3);
SU = std(mean(STORE,2),1,3);
LB = [-ones(prod(sz),1);-10*ones(sz(1),1);.01*ones(sz(1),1)];
UB = [ones(prod(sz),1);10*ones(sz(1),1);10*ones(sz(1),1)];
LB = [-ones(prod(sz),1);UU-10*SU;SU*.0001];
UB = [ones(prod(sz),1);UU+10*SU;SU/2];
xo = [rand(1,prod(sz)) 0*ones(1,sz(1)) 1*ones(1,sz(1))];
xTrain = particleswarm(func,numel(xo),LB,UB,options);
%%
global BD
BD = [];
sz = [3 3];
func = @(vec)countOP4(STORE,vec,NN,TARN,[],sz);
OF = @(A,B,C)outputPLT(func,A,B,C);
options = optimoptions('fmincon','Display','iter','PlotFcn',OF,'UseParallel',true);
xo = [rand(1,prod(sz)) 0*ones(1,sz(1)) 1*ones(1,sz(1))];
LB = [-ones(prod(sz),1);-10*ones(sz(1),1);.01*ones(sz(1),1)];
UB = [ones(prod(sz),1);10*ones(sz(1),1);10*ones(sz(1),1)];
xTrain = fmincon(func,xo,[],[],[],[],LB,UB,[],options);
%% try count
tmp = S1(:,:,e);
for e = 1:10
    %{
    tmpN = randi([15,30],1,1);
    tmp = generateImageP(MASK,R,tmpN);
    H1 = im2colF(double(tmp),[51 51],[ds ds]);
    H1 = sort(H1,1);
    Sc1 = PCA_REPROJ_T(H1,E1,U1);
    Sc1 = bsxfun(@times,Sc1,L1.^-1);
    %}
    Sc1 = STORE(:,:,e);
    tmpN = NN(e);
    [~,~,tmpCNT] = countOP4(Sc1,xTrain,25,[],[],sz);
    tmpCNT
    tmpN
end
%% next with real data - not yet
global BD
BD = [];
sz = [3 3];
func = @(vec)countOP4(STORE,vec,25,TARN,[],sz);
OF = @(A,B,C)outputPLT(func,A,B,C);
options = optimoptions('fmincon','Display','iter','PlotFcn',OF);
xo = [rand(1,prod(sz)) 0*ones(1,sz(1)) 1*ones(1,sz(1))];
LB = [-ones(prod(sz),1);-10*ones(sz(1),1);.01*ones(sz(1),1)];
UB = [ones(prod(sz),1);10*ones(sz(1),1);10*ones(sz(1),1)];
xTrain = fmincon(func,xo,[],[],[],[],LB,UB,[],options);
%%
[CNT,IMG] = countOP3(SCORE,xTrain,25,sz);
for k = 1:size(IMG,1)
    IMGV(:,:,k) = col2im(IMG(k,:),[51 51],size(tmp));
end
%%
func = @(vec)countOP2(STORE,vec,25);
T = reshape(STORE,[size(STORE,1) size(STORE,2)*size(STORE,3)]);
kidx = kmeans(T',2);
lda = myLDA(T',kidx);
scores = lda'*T;
for uTrain = 1:2
    U(uTrain) = mean(scores(kidx==uTrain));
    STD(uTrain) = std(scores(kidx==uTrain));
end
vec = [[U STD] lda'];
options = optimoptions('particleswarm','Display','iter');
xTrain = fminsearch(func,vec);
%% beta for counts on vec1 only
beta = vec1\ones(size(vec1,1),1);
%% beta for counts on vec2 and vec1
VEC = [vec1 vec2];
beta = VEC\ones(size(VEC,1),1);

%%
net = feedforwardnet(1);
view(net)

%%
[wS wC wU wE wL wERR wLAM] = PCA_FIT_FULL_T(W,3);
GMModel = fitgmdist(wC',3);
%%
wSZ = [51 51];
for e = 1%1:1:size(S,3)
    tic
    %tmp = S(:,:,e);
    [tmp] = imcrop(S(:,:,e),BOX);
    %mS(:,:,e) = imresize(S(:,:,e),1);
    H = im2colF(double(tmp)/255,wSZ,[1 1]);
    H1 = sort(H,1);
    
    [tmp] = imcrop(sD(:,:,e),BOX);
    %mS(:,:,e) = imresize(S(:,:,e),1);
    H = im2colF(double(tmp)/255,wSZ,[1 1]);
    H = sort(H,1);
    
    
    H = PCA_REPROJ_T([H;H1],wE,wU);
    [IDX,nlogl,P] = GMModel.cluster(H');
    IDX = col2im(IDX,wSZ,size(tmp));
    PK = col2im(P(:,3),wSZ,size(tmp));
end
%%
close all
for f = 1:size(S,3)
    imshow(S(:,:,f),[]);
    title(num2str(f))
    drawnow
end
%%
S = double(S)/255;
%%
close all
MASK = S(:,:,1) > graythresh(S(:,:,1))*.8;
imshow(MASK,[]);
drawnow
figure;
imshow(S(:,:,1),[]);
fidx = find(MASK);
%%
close all
imshow(mean(S,3),[]);
%% 
S0 = S(:,:,45:end);
uS = mean(S0,3);
uS = imfilter(uS,fspecial('gaussian',[51 51],9));
S0 = bsxfun(@minus,S0,uS);
imshow(uS,[]);
%%
close all
[n1,n2] = ndgrid(100:20:1000,150:20:1100);
pointList = [n1(:) n2(:)];
[pointListL pointListE SE SL] = wholeTrack_mod0(S0,pointList,1,[]);
%%
close all
for e = 1:size(S0,3)
    imshow(S0(:,:,e),[]);
    hold on
    plot(pointListL(:,2,e),pointListL(:,1,e),'r.')
    hold off
    drawnow
end
%%
pointListEM = pointListE;
pointListEM(:,:,1) = 0;
%%
clear X
[X(:,2) X(:,1) v] = impixel(S0(:,:,1),[]);
for t = 1:size(pointListE,3)
    F1{t} = griddedInterpolant(n1,n2,reshape(pointListEM(:,1,t),size(n1)));
    F2{t} = griddedInterpolant(n1,n2,reshape(pointListEM(:,2,t),size(n1)));
end
close all
clear X
MASK = S(:,:,1) > graythresh(S(:,:,45));
[X(:,1) X(:,2)] = find(MASK);
X = [n1(:) n2(:)];
for t = 1:size(pointListE,3)
    imshow(S0(:,:,t),[]);
    hold on
    plot(X(:,2,t),X(:,1,t),'r.')
    hold off
    drawnow
    d1 = F1{t}(X(:,1),X(:,2));
    d2 = F2{t}(X(:,1),X(:,2));
    X(:,1,t+1) = X(:,1,t) + d1;
    X(:,2,t+1) = X(:,2,t) + d2;
end
%%
close all
SKIP = 15;
[m1,m2] = ndgrid(100:SKIP:1000,150:SKIP:1100);
X = [m1(:) m2(:)];
d = 0;
for t = 1:size(pointListE,3)
    F1{t} = griddedInterpolant(n1,n2,reshape(pointListEM(:,1,t),size(n1)),'spline');
    F2{t} = griddedInterpolant(n1,n2,reshape(pointListEM(:,2,t),size(n1)),'spline');
end

SKIP = 1;
[xm1,xm2] = ndgrid(100:SKIP:1000,150:SKIP:1100);
XM = [xm1(:) xm2(:)];
d = 0;
for t = 1:size(pointListE,3)
    F1{t} = griddedInterpolant(n1,n2,reshape(pointListEM(:,1,t),size(n1)),'spline');
    F2{t} = griddedInterpolant(n1,n2,reshape(pointListEM(:,2,t),size(n1)),'spline');
end


PM = zeros(size(xm1));
oX = X;
oXM = XM;
V = [];
VELT = 1.5;



h1 = figure;
h2 = figure;
dataI = zeros(size(xm2,1),size(xm2,2),size(pointListE,3),3);
for t = 1:size(pointListE,3)
    Y = ba_interp2(S(:,:,t+44),xm2,xm1);
    YW = S(:,:,t+45);
   
    %CV = im2col(YW,[bSZ bSZ],'sliding');
    %[q1,q2] = ndgrid(-20:1:20,-20:1:20);
    
    d1 = F1{t}(X(:,1),X(:,2));
    d2 = F2{t}(X(:,1),X(:,2));
    VEL = sum([d1 d2].*[d1 d2],2).^.5;
    V = [V;VEL];
    %d1 = reshape(d1,size(m1));
    %d2 = reshape(d2,size(m1));
    DELTA = oX - X;
    X = X + [d1 d2];
   
    N = sum(DELTA.*DELTA,2).^.5;
    MAG = sum(DELTA.*DELTA,2);
    N(N == 0) = 1;
    DELTA  = bsxfun(@times,DELTA,N.^-1);
    X = X + .1*bsxfun(@times,DELTA,N);
    
    toRED = VEL > VELT;
    
    
    
    xd1 = F1{t}(XM(:,1),XM(:,2));
    xd2 = F2{t}(XM(:,1),XM(:,2));
    VELM = sum([xd1 xd2].*[xd1 xd2],2).^.5;
    DELTAM = oXM - XM;
    XT = XM;
    XM = XM + [xd1 xd2];
   
    NM = sum(DELTAM.*DELTAM,2).^.5;
    MAG = sum(DELTAM.*DELTAM,2);
    
    NM(NM == 0) = 1;
    DELTAM  = bsxfun(@times,DELTAM,NM.^-1);
    XM = XM + .1*bsxfun(@times,DELTAM,NM);
    
    dT = XM - XT;
    
    PM_TMP = zeros(size(YW));
    
    XM = round(XM);
    XM((XM(:,1) <= 0),1) = 1;
    XM((XM(:,1) > size(YW,1)),1) = size(YW,1);
    XM((XM(:,2) <= 0),2) = 1;
    XM((XM(:,2) > size(YW,2)),2) = size(YW,2);
    
    
    IDX = sub2ind(size(PM_TMP),XM(:,1),XM(:,2));
    
    %diskV1 = imfilter(reshape(xd1,size(xm1)),fspecial('disk',bSZ));
    %diskV2 = imfilter(reshape(xd2,size(xm1)),fspecial('disk',bSZ));
    %diskVT = imfilter(reshape(VELM,size(xm1)),fspecial('disk',bSZ));
    %toQ = diskVT > 5;
    
    
    
    %VELM = imfilter(reshape(VELM,size(xm1)),fspecial('gaussian',[31 31],7));
    %toPROB = VELM > 5.5;
    %toPROB = toPROB(:);
    %PM_TMP = PM;
    %PM_TMP(find(toPROB)) = 1;
    %PM_TMP(IDX(toPROB)) = 1;
    %PM_TMP = imfill(PM_TMP,'holes');
    %dB = bwboundaries(logical(PM_TMP));
    
    
    
    velR = X(toRED,:);
    dd = divergence(reshape(dT(:,2),size(xm1)),reshape(dT(:,1),size(xm2)));
    %dd = divergence(reshape(XM(:,2),size(xm1)),reshape(XM(:,1),size(xm2)));
    dd = imfilter(dd,fspecial('disk',bSZ),'replicate');
    %d = d + dd;
    d = dd;
    %Y = cat(3,bindVec(d),Y,Y);
    %imshow(Y,[]);
    
    tmp = ba_interp2(S(:,:,t+45),XM(:,2),XM(:,1));
    dataI(:,:,t,1) = reshape(tmp,size(xm1));
    dataI(:,:,t,2) = dd;
    dataI(:,:,t,3) = reshape(NM,size(xm1));
    %{
    tmp = permute(dataI(:,:,t,:),[1 2 4 3]);
    for k = 1:3
        tmp(:,:,k) = bindVec(tmp(:,:,k));
    end
    imshow(tmp,[]);
    %}
    
    
    %{
    figure(h1)
    imshow(d,[0 5]);
    hold on
    plot(X(:,2)-150,X(:,1)-100,'b.')
    plot(R(:,2)-150,R(:,1)-100,'r.')
    hold off
    %}
    
    %quiver(xm2(toQ),xm1(toQ),diskV2(toQ),diskV1(toQ),10)
    %   hold off
    
    %{
    figure(h2);
    O = flattenMaskOverlay(YW,logical(PM_TMP));
    imshow(O,[]);
    hold on
    %for c = 1:numel(dB)
    %    plot(dB{c}(:,2),dB{c}(:,1),'r')
    %end
    hold off
    %}
    %waitforbuttonpress
    %pause(.5)
    %hello=1
    
    t
    %drawnow
end
%%
