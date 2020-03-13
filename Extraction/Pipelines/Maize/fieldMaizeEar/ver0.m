FilePath = '/mnt/spaldingdata/nate/octerineDataStore/maizeData/fieldEarData/ear_images';
FileList = {};
FileExt = {'jpg'};
FileList = fdig(FilePath,FileList,FileExt,1);
%%
for e = 1:numel(FileList)
    I = imread(FileList{e});
    imshow(I,[]);
    drawnow
end
%%
close all
areaThresh = 4000;
oPath = '/mnt/spaldingdata/nate/octerineDataStore/maizeData/fieldEarData/ear_images/return2/';
mkdir(oPath);

for e = 1:numel(FileList)
    [pth,nm,ext] = fileparts(FileList{e});
    
    
    border_oFile = [oPath nm '_border.csv'];
    central_oFile = [oPath nm '_central.csv'];
    image_oFile = [oPath nm '_return.jpg'];
    I = imread(FileList{e});
    
    %{
    if rand(1) > .5
        I = imrotate(I,180);
    end
    %}
    LAB = rgb2lab(I);
    
    sig0 = LAB(:,:,1);
    
    sig = LAB(:,:,3);
    sig = bindVec(sig);
    th = graythresh(sig);
    
    sig0 = bindVec(sig0);
    th0 = graythresh(sig0);
    whiteBoard = sig0 > th0;
    whiteBoard = imclose(whiteBoard,strel('disk',21,0));
    whiteBoard = imfill(whiteBoard,'holes');
    whiteBoard = bwlarge(whiteBoard);
    MASK = sig > th & whiteBoard;
    
    
    
    MASK = bwareaopen(MASK,areaThresh);
    
    R = regionprops(MASK,'BoundingBox','MajorAxis','MinorAxis','Area','Centroid');
    
    topV = [];
    for e = 1:numel(R)
        topV(e) = R(e).Centroid(2) > size(I,1)/2;
    end
    
    deltaTOP = [];
    deltaBOT = [];
    for e = 1:numel(R)
        if topV(e)
            deltaTOP = [deltaTOP R(e).Centroid(1)];
        else
            deltaBOT = [deltaBOT R(e).Centroid(1)];
        end
    end
    deltaTOP = sort(deltaTOP);
    deltaBOT = sort(deltaBOT);
    uT = mean(diff(deltaTOP));
    uB = mean(diff(deltaBOT));
    
    if uT > uB
        topBorder = true;
    else
        topBorder = false;
    end
    
    
    
    out = flattenMaskOverlay(I,MASK);
    imshow(out,[]);
    hold on
    
    borderName = {};
    centralName = {};
    
    borderArea = [];
    centralArea = [];
    
    borderLength = [];
    centralLength = [];
    
    borderWidth = [];
    centralWidth = [];
    
    
    for e = 1:numel(R)
        rectangle('Position',R(e).BoundingBox,'EdgeColor','g');
        L1(1,:) = R(e).BoundingBox(1:2);
        L1(2,:) = R(e).BoundingBox(1:2) + [R(e).BoundingBox(3) 0];
        L1(:,2) = L1(:,2) - 40;
        
        L2(1,:) = R(e).BoundingBox(1:2);
        L2(2,:) = R(e).BoundingBox(1:2) + [0 R(e).BoundingBox(4)];
        L2(:,1) = L2(:,1) - 40;
        
        plot(L1(:,1),L1(:,2),'b') 
        plot(L2(:,1),L2(:,2),'b') 
        
        T1 = L1(1,:) - 50;
        text(T1(1),T1(2),num2str(round(R(e).MinorAxisLength)));
        
        T2 = L2(1,:) - 20;
        h = text(T2(1),T2(2),num2str(round(R(e).MajorAxisLength)));
        set(h,'Rotation',90+180);
        
        
        if topBorder     
            if topV(e)
                nCL = 'm';
                borderArea = [borderArea R(e).Area];
                borderLength = [borderLength R(e).MajorAxisLength];
                borderWidth = [borderWidth R(e).MinorAxisLength];
                borderName{end+1} = nm;
            else
                nCL = 'c';
                centralArea = [centralArea R(e).Area];
                centralLength = [centralLength R(e).MajorAxisLength];
                centralWidth = [centralWidth R(e).MinorAxisLength];
                centralName{end+1} = nm;
            end
        else
            if topV(e)
                nCL = 'c';
                centralArea = [centralArea R(e).Area];
                centralLength = [centralLength R(e).MajorAxisLength];
                centralWidth = [centralWidth R(e).MinorAxisLength];
                centralName{end+1} = nm;
            else
                nCL = 'm'; 
                borderArea = [borderArea R(e).Area];
                borderLength = [borderLength R(e).MajorAxisLength];
                borderWidth = [borderWidth R(e).MinorAxisLength];
                borderName{end+1} = nm;
            end
        end
        dX = [100 -100];
        text(T1(1)+dX(1),T1(2)+dX(2),num2str(e),'BackgroundColor',nCL);
    end
    hold off
    saveas(gca,image_oFile);
    
    
    borderTable = table(borderName',borderArea',borderLength',borderWidth','VariableNames',{'Name','Area','Length','Width'});
    centralTable = table(centralName',centralArea',centralLength',centralWidth','VariableNames',{'Name','Area','Length','Width'});
    
   
    writetable(borderTable,border_oFile);
    writetable(centralTable,central_oFile);
     
    drawnow
end
