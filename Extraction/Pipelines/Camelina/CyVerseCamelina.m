function [] = CyVerseCamelina(fileName,oPath)

    %     iter = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% Prepare the Mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CamelinaIm = imread(fileName);
    %       figure, imshow(CamelinaIm); impixelinfo;  
    % hard code for creating mask - comment added via nate miller
    % removing "small objects"
    Mask01 = CamelinaIm(:, :, 1) >= 30 & CamelinaIm(:, :, 2) >= 30 & CamelinaIm(:, :, 3) >= 25; 
    areaRemoved = bwareaopen(Mask01, 10000);
    %     figure, imshow(areaRemoved); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% label and measure each object %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % lable each object
    [L3, CamelinaNum3]=bwlabel(areaRemoved,4);
    %     figure, imshow(L, []); impixelinfo;
%     % Obtatin geometric features
    CamelinaStat = regionprops(L3,'Area','MinorAxisLength');  % 



%     areas = [CamelinaStat.Area]; 
%     [maxArea maxAreaIdx]= max(areas); 
%     sort(areas, 'descend');
    
    % Calculate pixel to actual size coefficient
    % quarter coin size: 24.26 mm; dime coin size: 17.91 mm
    %sizeConst = 17.91/CamelinaStat(maxAreaIdx).MajorAxisLength;
    % Remove small objects that are not Camelinas and reference coins

    
    for i = 1:CamelinaNum3
        if (CamelinaStat(i).Area > 10000 && CamelinaStat(i).Area < 200000)
            sizeConst = 17.91/CamelinaStat(i).MinorAxisLength;
            index = find(L3==i);
            L3(index) = 0;  
        end 
    end 
    
    Mask02 = abs(CamelinaIm(:, :, 1) - CamelinaIm(:, :, 2)) < 20; 
    %     figure, imshow(Mask01);  figure, imshow(Mask02);  
%     areaRemoved02 = bwareaopen(Mask02, 500);
%     figure, imshow(areaRemoved02); 

%     % Morphological operations
%     Mask_Close = bwmorph(Mask, 'close'); 
% %     figure, imshow(Mask_Close);

    Mask_Fill = bwmorph(Mask02, 'fill',5); 
%     figure, imshow(Mask_Fill, []);

%     SE = strel('disk', 7);
%     Mask_Open = imopen(Mask_Fill, SE);
%     figure, imshow(Mask_Open); impixelinfo; 
    
%     SE2 = strel('disk', 3);
%     Mask_Dilate = imdilate(~Mask_Fill, SE2);
%     figure, imshow(Mask_Dilate); impixelinfo; 

    SE2 = strel('disk', 1);
    Mask_Erode = imerode(~Mask_Fill, SE2);
    areaRemoved02 = bwareaopen(Mask_Erode, 500);
%     figure, imshow(areaRemoved02); 
%     figure, imshow(areaRemoved02); impixelinfo; 


%     figure, imshow(CamelinaIm);
%     hold on; 
%     IMhandle = imshow(Mask_Erode);
%     set(IMhandle, 'AlphaData', 0.3);
%     hold off;
%     impixelinfo;

%     IM = Mask_Erode; 
    
    [L, CamelinaNum]=bwlabel(Mask_Erode,4);
%     figure, imshow(L, []);

%     % Obtatin geometric features
    CamelinaStat = regionprops(L,'Area','MajorAxisLength');  % 

%     areas = [CamelinaStat.Area]; 
%     [maxArea maxAreaIdx]= max(areas); 
%     sort(areas, 'descend');
    
    % Calculate pixel to actual size coefficient
    % quarter coin size: 24.26 mm; dime coin size: 17.91 mm
    %sizeConst = 17.91/CamelinaStat(maxAreaIdx).MajorAxisLength;
    % Remove small objects that are not Camelinas and reference coins

    
    for i = 1:CamelinaNum
        
        if (CamelinaStat(i).Area > 10000 || CamelinaStat(i).Area < 500 || CamelinaStat(i).MajorAxisLength > 300)  % Area needs to be optimized. | CamelinaStat(i).Eccentricity > 0.6
             index = find(L==i);
             L(index) = 0;        
        end 
    end

    
    CamelinaMask = im2bw(L);
% %     clear CamelinaStat; 
%     figure, imshow(CamelinaMask,[]);
 
    IM = CamelinaMask;
 
    revIM = imcomplement(IM);
    %figure()
    %subplot (1,3,1)
    %imshow(revim,[])
    %title('imcomplement to watershed')

    D = bwdist(revIM);
    D2 = imcomplement(D);
    %subplot (1,3,2)
    %imshow(D2,[])
    %title('Distance image to watershed')
    % Suppress shallow minima
    D3 = imhmin(D2,1);
    L1 = watershed(D3);
    IM(L1 == 0) = 0;
%     figure, imshow(IM);
    %%%%%% Watershed operation ends %%%%%%
    
    % When analyzing a batch of image, uncomment the set() to stop showing
    % the resulted image
    figure, % set(gcf,'Visible','off'); % Do not show the image
    imshow(IM); impixelinfo;
    hold on; 
    IMhandle = imshow(CamelinaIm);
    set(IMhandle, 'AlphaData', 0.5);
    hold off;
    strTitle = strrep(fileName, '.JPG', ''); 
    strTitle = strrep(strTitle, '_', ' '); 
    title(strTitle);
    

 %%%%%%%%%%%%%%%%%%%%% Extract features %%%%%%%%%%%%%%%%%%%%%%%   
    [L2, CamelinaNum2]=bwlabel(CamelinaMask,4);
%     figure, imshow(L2, []);


%     % Obtatin geometric features
    CamelinaStat2 = regionprops(L2,'Area','MajorAxisLength','MinorAxisLength','Eccentricity');  
    
    CamelinaTbl = cell(CamelinaNum2, 5);
               
%   % Use file name as entry name
    CamelinaTbl(1:CamelinaNum2, 1) = {strrep(fileName, '.JPG', '')};
    % Number of Camelina in one image
    for m = 1 : CamelinaNum2
        CamelinaTbl(m,2) = {CamelinaStat2(m).Area*sizeConst*sizeConst}; 
        CamelinaTbl(m,3) = {CamelinaStat2(m).MajorAxisLength*sizeConst}; 
        CamelinaTbl(m,4) = {CamelinaStat2(m).MinorAxisLength*sizeConst}; 
        CamelinaTbl(m,5) = {CamelinaStat2(m).Eccentricity};        
    end 
    
    
    fileStr = strcat(oPath,strrep(fileName, '.jpg', '.fig')); 
%     saveas(gcf, char(fileStr));

    
    fileStr = strcat(oPath,strrep(fileName, '.jpg', '.mat')); 
%     save(char(fileStr),'CamelinaTbl');
    
%%%%%%%%%%%%%%%%%%%%%%%% Store info for each variety %%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    sortedCamlinaTbl = cell(CamelinaNum2, 6);
    sortedCamlinaTbl = sortrows(CamelinaTbl, -2); 
    % Prepared to average the middle 80% data of 'MajorAxisLength','MinorAxisLength','Eccentricity'; 
    num10Per = floor(CamelinaNum2*0.1); 
    num90Per = floor(CamelinaNum2*0.9); 

    summaryTbl(1,1) ={ fileName };
    summaryTbl(1,2) ={ CamelinaNum2 };
    summaryTbl(1,3) ={ mean(cell2mat(sortedCamlinaTbl(num10Per:num90Per, 2))) };
    summaryTbl(1,4) ={ mean(cell2mat(sortedCamlinaTbl(num10Per:num90Per, 3))) };
    summaryTbl(1,5) ={ mean(cell2mat(sortedCamlinaTbl(num10Per:num90Per, 4))) };
    summaryTbl(1,6) ={ mean(cell2mat(sortedCamlinaTbl(num10Per:num90Per, 5))) };
       

end

%{



    fileName = '/home/nate/11603.jpg';
    %fileName = '/home/nate/412_36_10.jpg';
    oPath = './output/';
    %CyVerseCamelina(fileName,oPath);
    measureArabidopsisSeed(fileName,oPath);



%}
