pvcInnerRadius = 120;
pvcThickness = 27;
pvcOuterRadius = pvcInnerRadius + pvcThickness;
mag = 4;
capRadius = 115;
capHeight = 20;
height = 300;
RAD = 2000;
PIPE_N = [RAD,1,300];
RIM_N = [RAD,mag*pvcThickness,1];
CAPTOP_N = [RAD,mag*capRadius,1];
CAPRIM_N = [RAD,1,capHeight];
%% inner cyliner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x1,x2,x3] = ndgrid(linspace(-pi,pi,PIPE_N(1)),...
                    linspace(pvcInnerRadius,pvcInnerRadius,PIPE_N(2)),...
                    linspace(0,height,PIPE_N(3)));
sur1_1 = x2.*cos(x1);
sur1_2 = x2.*sin(x1);
sur1_3 = x3;
sur1_1 = squeeze(sur1_1);
sur1_2 = squeeze(sur1_2);
sur1_3 = squeeze(sur1_3);
surface(sur1_1,sur1_2,sur1_3);hold on
%% outer cyliner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x1,x2,x3] = ndgrid(linspace(-pi,pi,PIPE_N(1)),...
                    linspace(pvcOuterRadius,pvcOuterRadius,PIPE_N(2)),...
                    linspace(0,height,PIPE_N(3)));
sur2_1 = x2.*cos(x1);
sur2_2 = x2.*sin(x1);
sur2_3 = x3;
sur2_1 = squeeze(sur2_1);
sur2_2 = squeeze(sur2_2);
sur2_3 = squeeze(sur2_3);
load clown
C = imresize(flipud(X),size(sur2_1));

surface(sur2_1,sur2_2,sur2_3,C,'FaceColor','texturemap',...
    'EdgeColor','none',...
    'CDataMapping','direct');
%% rim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x1,x2,x3] = ndgrid(linspace(-pi,pi,RIM_N(1)),...
                    linspace(pvcInnerRadius,pvcOuterRadius,RIM_N(2)),...
                    linspace(0,0,RIM_N(3)));
sur3_1 = x2.*cos(x1);
sur3_2 = x2.*sin(x1);
sur3_3 = x3;
sur3_1 = squeeze(sur3_1);
sur3_2 = squeeze(sur3_2);
sur3_3 = squeeze(sur3_3);
surface(sur3_1,sur3_2,sur3_3);
%% cap top
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x1,x2,x3] = ndgrid(linspace(-pi,pi,CAPTOP_N(1)),...
                    linspace(0,capRadius,CAPTOP_N(2)),...
                    linspace(0,0,CAPTOP_N(3)));
sur4_1 = x2.*cos(x1);
sur4_2 = x2.*sin(x1);
sur4_3 = x3;
sur4_1 = squeeze(sur4_1);
sur4_2 = squeeze(sur4_2);
sur4_3 = squeeze(sur4_3);
surface(sur4_1,sur4_2,sur4_3);hold on
%% cap rim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x1,x2,x3] = ndgrid(linspace(-pi,pi,CAPRIM_N(1)),...
                    linspace(capRadius,capRadius,CAPRIM_N(2)),...
                    linspace(0,capHeight,CAPRIM_N(3)));
sur5_1 = x2.*cos(x1);
sur5_2 = x2.*sin(x1);
sur5_3 = x3;
sur5_1 = squeeze(sur5_1);
sur5_2 = squeeze(sur5_2);
sur5_3 = squeeze(sur5_3);
surface(sur5_1,sur5_2,sur5_3);
%% make movie version 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
capImage = imread('/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/Maize/maizeEmergence/modelMethod/modelMean.tif');
capImage = double(capImage)/255;


FilePath = '/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/Maize/maizeEmergence/modelMethod/imodels/';
FileList = {};
FileExt = {'tif'};
tic
FileList = gdig(FilePath,FileList,FileExt,1);
capImage = double(imread(FileList{1}))/255;
imshow(capImage,[]);
hold on
%plot(sur4_1(:)+150,sur4_2(:)+150,'.')
%plot(sur3_1(:)+150,sur3_2(:)+150,'.')
capTexture = squeeze(ba_interp2(capImage,sur4_1+150,sur4_2+150));
rimTexture = ba_interp2(capImage,sur3_1+150,sur3_2+150);
%% make surface with texture
close all
surface(sur4_1,sur4_2,sur4_3,capTexture,'FaceColor','texturemap',...
    'EdgeColor','none',...
    'CDataMapping','direct');
hold on
surface(sur3_1,sur3_2,sur3_3,rimTexture,'FaceColor','texturemap',...
    'EdgeColor','none',...
    'CDataMapping','direct');
surface(sur2_1,sur2_2,sur2_3,'EdgeColor','none');
surface(sur1_1,sur1_2,sur1_3,'EdgeColor','none');
surface(sur5_1,sur5_2,sur5_3);
view([0 -90]);





