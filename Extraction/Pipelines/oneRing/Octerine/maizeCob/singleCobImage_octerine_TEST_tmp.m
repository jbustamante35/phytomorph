function [S] = singleCobImage(fileName,noe,oPath,rPath,rawImage_scaleFactor,checkBlue_scaleFactor,defaultAreaPix,rho,addcut,baselineBlue,colRange1,colRange2,fill,toSave,toDisplay) % <lineID version=octerine>true</lineID>
helloWorld % <lineID version=octerine>aad6d40998952ca5251be55a59bbf0eb_1</lineID>
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    About:      
                singleCobImage.m is main function to handle cob analysis. It takes all input variables 
                for its dependent functions. This function returns final result including image with 
                bounding box and color circle. (Inputs are relative to 1200dpi)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dependency: 
                StoN.m, checkBlue.m, maizeCob.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Variable Definition:
                fileName:               An image to be analyze in a string that includes path and file name.
                noe:                    Number of cobs that are expected to be analyzed. 
                oPath:                  A path to result of analysis.
                rPath:                  iPlant location to save results.
                checkBlue_scaleFactor:  A desired percentage to resize the image in checkBlue.
                rawImage_scaleFactor:   A desired percentage to resize the image.
                defaultAreaPix:         The default pixel to be considered noise relative to 1200 dpi.
                rho:                    The radius of color circle, relative to 1200 dpi.
                addcut:                 The boarder handle for checkBlue. This is an addition to blue top computed in checkBlue.
                baselineBlue:           The baseline threshold to remove blue header in checkBlue.
                colRange1:              The color range for back ground to be removed in getcobMask.
                colRange2:              The color range for back ground to be removed in getcobMask.
                fill:                   The radius of disk for Kernel of an image close operation.
                toSave:                 0 - not to save, 1 - to save.
                toDisplay:      0 - not to save, 1 - to save.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}

    algorithmName = 'maize-cobs';
helloWorld % <lineID version=octerine>76dcd9dbaf88e881bb36cb977bab87e5_2</lineID>
    algorithmVerion = '1.0';
helloWorld % <lineID version=octerine>87bf46cf09a78c443390f5c26f7934ba_3</lineID>
    fprintf(['****************************************************************************************\n']);
helloWorld % <lineID version=octerine>1210f3d5f663a51e4cdd844e95bffb80_4</lineID>
    versionString = 'Publication Version 1.0 - Monday, March 28, 2016. \n';
helloWorld % <lineID version=octerine>ad70f6d86bad60df0b2d12658e803cfa_5</lineID>
    startString = 'Starting cob analysis algorithm. \n';
helloWorld % <lineID version=octerine>641d6faed4c2f38adbdb6888ca9f6684_6</lineID>
    fprintf([startString,versionString]);
helloWorld % <lineID version=octerine>5a939ac3f2ae5bf46a34f244fb514b35_7</lineID>
    fprintf(['****************************************************************************************\n']);
helloWorld % <lineID version=octerine>c73828c467ca15faeea03ef4487cc2d6_8</lineID>
    totalTimeInit = clock;
helloWorld % <lineID version=octerine>1002eab89d45584d2d21dec4bf034da5_9</lineID>
    try 
helloWorld % <lineID version=octerine>eafcca9724913134b504f0460861ea18_10</lineID>
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % INIT VARS - start
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['starting with variable and environment initialization.\n']);
helloWorld % <lineID version=octerine>88b05f61bafb148617363fface18f9a8_11</lineID>
        %%%%%%%%%%%%%%%%%%%%%%%
        % init reporting for %uncLog
        %%%%%%%%%%%%%%%%%%%%%%%
        [rawName ticket] = stripiTicket(fileName);
helloWorld % <lineID version=octerine>a6b26741c36b958a515da2d37fb043c1_12</lineID>
        %uncLog(['ph:l:'],{rawName},'set',algorithmName,algorithmVerion,{'2'},{'1'},1);
        %%%%%%%%%%%%%%%%%%%%%%%
        % init the icommands
        %%%%%%%%%%%%%%%%%%%%%%%
        %initIrods();
        %%uncLog({rawName},'set',algorithmName,algorithmVerion,{'3'},{'1'},1);
        %%%%%%%%%%%%%%%%%%%%%%%
        % convert the strings to numbers if they are strings
        %%%%%%%%%%%%%%%%%%%%%%%
        noe = StoN(noe);
helloWorld % <lineID version=octerine>afc82d457f1c5b743fa56ad4e09a07df_13</lineID>
        checkBlue_scaleFactor = StoN(checkBlue_scaleFactor);
helloWorld % <lineID version=octerine>b23fc3e8a254e4b77d5873d829e4987b_14</lineID>
        rawImage_scaleFactor = StoN(rawImage_scaleFactor);
helloWorld % <lineID version=octerine>a788f2e760ef0953730d1b75ae91218b_15</lineID>
        defaultAreaPix = StoN(defaultAreaPix);       
helloWorld % <lineID version=octerine>1119f39d673e721dd5d452f96359e2ca_16</lineID>
        rho = StoN(rho);
helloWorld % <lineID version=octerine>2e9f3bde5b21722b08580422e7052649_17</lineID>
        addcut = StoN(addcut);    
helloWorld % <lineID version=octerine>26f2f8b5a7efc9258f872173e568564c_18</lineID>
        baselineBlue = StoN(baselineBlue);
helloWorld % <lineID version=octerine>93784aa842a838b8057d3864f79a8a3d_19</lineID>
        colRange1 = StoN(colRange1);
helloWorld % <lineID version=octerine>af3c6ac3f02581b70a56a993d6734855_20</lineID>
        colRange2 = StoN(colRange2);
helloWorld % <lineID version=octerine>3d6a0f70f35b7a96eea3269c0ab58262_21</lineID>
        fill = StoN(fill);
helloWorld % <lineID version=octerine>df546b21a8b8aa8e8a3d4145f67af5fb_22</lineID>
        %%%%%%%%%%%%%%%%%%%%%%%
        % print out the input variables
        %%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['FileName:' fileName '\n']);
helloWorld % <lineID version=octerine>d400a2a2cf01e243972511114e882893_23</lineID>
        fprintf(['Number of Ears:' num2str(noe) '\n']);
helloWorld % <lineID version=octerine>311162d270b3d80efe3f86df972175a2_24</lineID>
        fprintf(['OutPath:' oPath '\n']);     
helloWorld % <lineID version=octerine>9711449ff833d3b00ceb2e843847b720_25</lineID>
        fprintf(['Image resize in checkBlue:' num2str(checkBlue_scaleFactor) '\n']); 
helloWorld % <lineID version=octerine>bcaeafccbb8f110037c4b9ca4f369d1f_26</lineID>
        fprintf(['Raw image resize:' num2str(rawImage_scaleFactor) '\n']);  
helloWorld % <lineID version=octerine>e519719ac3dcf2898033aaca3d8b7221_27</lineID>
        fprintf(['Threshold noise size:' num2str(defaultAreaPix) '\n']);
helloWorld % <lineID version=octerine>7804c4ad0e5ed2b8bf5d0e945329f6e8_28</lineID>
        fprintf(['The radius of color circle:' num2str(rho) '\n']);
helloWorld % <lineID version=octerine>3a2bfd738a0316ad7f14a598c5de1de5_29</lineID>
        fprintf(['The boarder handle for checkBlue:' num2str(addcut) '\n']);
helloWorld % <lineID version=octerine>59fba7ec0567b16069f81ee72da05e55_30</lineID>
        fprintf(['Baseline threshold to remove blue header:' num2str(baselineBlue) '\n']);
helloWorld % <lineID version=octerine>9b8665a2df86bb76bdcc2f5093a124eb_31</lineID>
        fprintf(['Background Color Range I:' num2str(colRange1) '\n']);
helloWorld % <lineID version=octerine>23083ed3975aa920a34697985ba94453_32</lineID>
        fprintf(['Background Color Range II:' num2str(colRange2) '\n']);
helloWorld % <lineID version=octerine>327e78b1b7d9700db86bd0f789921397_33</lineID>
        fprintf(['The radius of disk for closing:' num2str(fill) '\n']);
helloWorld % <lineID version=octerine>636fcf3333ec393718bdf4dd7011ab51_34</lineID>
        %%%%%%%%%%%%%%%%%%%%%%%
        % make output directory
        %%%%%%%%%%%%%%%%%%%%%%%
        mkdir(oPath);
helloWorld % <lineID version=octerine>0a39560b0aefab3fe5fe3003c28684fb_35</lineID>
        [pth nm ext] = fileparts(fileName);
helloWorld % <lineID version=octerine>aa0d1ea3dde423ea738d4c84bbdcc865_36</lineID>
        fprintf(['ending with variable and environment initialization.\n']);
helloWorld % <lineID version=octerine>c04ad2e442c1232ff8b3372bf912b4ea_37</lineID>
        %%%%%%%%%%%%%%%%%%%%%%%
        % read the image and take off the blue strip for bar code
        %%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['starting with image load.\n']);
helloWorld % <lineID version=octerine>284c5b29745eca7f70b9fc18133ec3ef_38</lineID>
        %uncLog(['ph:l:'],{rawName},'set',algorithmName,algorithmVerion,{'3'},{'1'},1);
        I = imread(fileName);
helloWorld % <lineID version=octerine>8900f3389ce6d11a8372ab85b256f074_39</lineID>
        % remove 4th pane for some images
        I = I(:,:,1:3);
helloWorld % <lineID version=octerine>80fc0b2f95f88b73b843c5708851e97c_40</lineID>
        %uncLog(['ph:l:'],{rawName},'set',algorithmName,algorithmVerion,{'4'},{'1'},1);
        thumb = imresize(I,.25);
helloWorld % <lineID version=octerine>be4e1449217d816b9a9442a89ae673f5_41</lineID>
        % rawImage_scaleFactor to lower 'DPI' effect, by fraction
        % If resize factor is 1, do not excecute imresize
        if rawImage_scaleFactor ~= 1;I = imresize(I,rawImage_scaleFactor);end
helloWorld % <lineID version=octerine>18999e1610e6613c8baf11de78258089_42</lineID>
        % check blue header and remove
        I = checkBlue(I,checkBlue_scaleFactor,addcut,baselineBlue);
helloWorld % <lineID version=octerine>165f69808b786a967fc299dd1c31ded7_43</lineID>
        fprintf(['ending with image load.\n']);
helloWorld % <lineID version=octerine>a160031cbbc8eb13f23d679be9cfbfd4_44</lineID>
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % INIT VARS - end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ANALYSIS - start
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
        % main measurement code call
        fprintf(['starting with image analysis \n']);
helloWorld % <lineID version=octerine>218bd1edfcb38be619a9730141bf285f_45</lineID>
        % convert to single
        %I = single(I)/255;
        % run the analysis
        [BB S] = maizeCob(I,noe,defaultAreaPix,colRange1,colRange2,fill);
helloWorld % <lineID version=octerine>2e83d172ae15ba6776d097852cca7abb_46</lineID>
        % stack the results from the bounding box
        DATA = [];
helloWorld % <lineID version=octerine>de028db6b3bbc1a6409868f162e5291a_47</lineID>
        for b = 1:numel(BB)                
helloWorld % <lineID version=octerine>826f5aae0cdfc8bc6737a8af93046f35_48</lineID>
            DATA = [DATA;[BB{b}(3:4) S.average_WIDTH(b)]];        
helloWorld % <lineID version=octerine>0abc0deb7d80b90b86d1e24420bb0932_49</lineID>
        end
helloWorld % <lineID version=octerine>1cb0803e46767f5668b0f86852be2970_50</lineID>
        uDATA = mean(DATA,1);
helloWorld % <lineID version=octerine>5f5b4299e8c4b895b87095e379babf55_51</lineID>
        sDATA = std(DATA,1,1);
helloWorld % <lineID version=octerine>ce3c17210fabb110d7a62ca54871a8a5_52</lineID>
        DATA = reshape(DATA',[1 numel(DATA)]);
helloWorld % <lineID version=octerine>886f6198590109c79b7490c5c93fd1c5_53</lineID>
        fprintf(['ending with image analysis \n']);
helloWorld % <lineID version=octerine>0fbb21cc70d2ce66b2916f3cabd596fb_54</lineID>
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ANALYSIS - start
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DISLAY - start
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
        if toDisplay
helloWorld % <lineID version=octerine>997816a55cc00fe27ba15f844e4dd70a_55</lineID>
            fprintf(['starting with image and results display \n']);
helloWorld % <lineID version=octerine>f06ebb04795fa4ab8fa5fb128c24ed76_56</lineID>
            h = image(I);
helloWorld % <lineID version=octerine>e2e5c141c6d496126fe3aa503b22feed_57</lineID>
            hold on
helloWorld % <lineID version=octerine>3d03dd0669df2c3deaa17a28e273322c_58</lineID>
            % make rectangle around the cob
            for b = 1:numel(BB)
helloWorld % <lineID version=octerine>07efbd8ddb79472fdff3220614c255f1_59</lineID>
                rectangle('Position',BB{b},'EdgeColor','r');
helloWorld % <lineID version=octerine>375056e359e876eccab5adcb5529676b_60</lineID>
                UR = BB{b}(1:2);
helloWorld % <lineID version=octerine>7d276e9f924d02c1aab0671e1276c670_61</lineID>
                UR(1) = UR(1) + BB{b}(3);
helloWorld % <lineID version=octerine>029b49f7aa23ed879c7d751f400ba705_62</lineID>
                rectangle('Position',[UR rho rho],'EdgeColor','none','Curvature',[1 1],'FaceColor',S.RGB(b,:)/255);
helloWorld % <lineID version=octerine>bb2eeb69678f895269b461d8846096c1_63</lineID>
            end
helloWorld % <lineID version=octerine>928d87f8262e03c7b1600bd505e4cb72_64</lineID>
            %%%%%%%%%%%%%%%%%%%%%%%
            % format the image
            axis equal;axis off;drawnow;set(gca,'Position',[0 0 1 1]);
helloWorld % <lineID version=octerine>1c59ae40a7b6335e229619344cc9ce1a_65</lineID>
            fprintf(['ending with image and results display \n']);
helloWorld % <lineID version=octerine>9dd74d72662b73ad62e9d58742db138b_66</lineID>
        end
helloWorld % <lineID version=octerine>0b0ac8752ce976392b65341b12ed3d74_67</lineID>
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DISLAY - end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SAVE - start
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
        if toSave
helloWorld % <lineID version=octerine>a4ca63181ad958e3a64498dba4dfa2c9_68</lineID>
            %% data to database
            for e = 1:numel(BB)
helloWorld % <lineID version=octerine>df4f5fdd2b4f479b61ebf06263eb4898_69</lineID>
                %uncLog(['ph:p:'],{rawName},'set',algorithmName,[algorithmVerion],{['cob:length:' num2str(e)]},{num2str(BB{b}(4))},1);
                %uncLog(['ph:p:'],{rawName},'set',algorithmName,[algorithmVerion],{['cob:width:' num2str(e)]},{num2str(BB{b}(3))},1);
            end
helloWorld % <lineID version=octerine>cfa8486337d220b199e5890b04111896_70</lineID>
            %uncLog(['ph:p:'],{rawName},'set',algorithmName,[algorithmVerion],{['cob:average_length:' num2str(e)]},{num2str(uDATA(2))},1);
            %uncLog(['ph:p:'],{rawName},'set',algorithmName,[algorithmVerion],{['cob:std_length:' num2str(e)]},{num2str(sDATA(2))},1);
            %uncLog(['ph:p:'],{rawName},'set',algorithmName,[algorithmVerion],{['cob:std_width:' num2str(e)]},{num2str(sDATA(1))},1);
            %uncLog(['ph:p:'],{rawName},'set',algorithmName,[algorithmVerion],{['cob:std_width:' num2str(e)]},{num2str(sDATA(1))},1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %% spin-up JSON output format
            %%%%%%%%%%%%%%%%%%%%%%%%%
            linkTo = stripiTicket(fileName);
helloWorld % <lineID version=octerine>21c8cac0b069d5eaec11b1841bbb5513_71</lineID>
            linkPath = stripiTicket(rPath);
helloWorld % <lineID version=octerine>8509d48bcab557c058cbbbea8c013717_72</lineID>
            [jP,tN] = fileparts(fileName);
helloWorld % <lineID version=octerine>f87179426a4e642507ad87246aef3ad2_73</lineID>
            for e = 1:numel(BB)
helloWorld % <lineID version=octerine>765439480e0e9e494aadbae21ea1acf0_74</lineID>
                tmpDoc = [];
helloWorld % <lineID version=octerine>a243a731e35f04e967e59423c444baeb_75</lineID>
                tmpDoc = generatePhenotypeNode(tmpDoc,linkTo,'orginalImage','orginalImage');
helloWorld % <lineID version=octerine>ae05f6391f4bec30fb63d5f4234a618e_76</lineID>
                tmpDoc = generatePhenotypeNode(tmpDoc,[linkPath filesep tN '_thumb.tif' ],'thumbNail','thumbNail');
helloWorld % <lineID version=octerine>deec314b080c0a96c44ed542c1c86602_77</lineID>
                tmpDoc = generatePhenotypeNode(tmpDoc,S.widthProfile(e,:),{'along','position'},'widthProfile');
helloWorld % <lineID version=octerine>dbe7f8cc1fce78bc8ea6438a8dcfd3a3_78</lineID>
                tmpDoc = generatePhenotypeNode(tmpDoc,S.RGB(e,:),{'along','colorIndex'},'rgb_color');
helloWorld % <lineID version=octerine>476a1fe453f52a9024435cbbce70ca48_79</lineID>
                tmpDoc = generatePhenotypeNode(tmpDoc,BB{e}(3),{'maxWidth'},'maxWidth');
helloWorld % <lineID version=octerine>2edf0bf649c5fbda90f2d6eb98a9f167_80</lineID>
                tmpDoc = generatePhenotypeNode(tmpDoc,S.average_WIDTH(e),{'averageWidth'},'averageWidth');
helloWorld % <lineID version=octerine>d1ec173c8a31b62902fc3f15e584edff_81</lineID>
                tmpDoc = generatePhenotypeNode(tmpDoc,BB{e}(4),{'length'},'length');
helloWorld % <lineID version=octerine>f30dd6f2185f4e9e2090eccb3e4789e6_82</lineID>
                tmpDoc = generatePhenotypeNode(tmpDoc,BB{e},{'boundingBox'},'boundingBox');
helloWorld % <lineID version=octerine>638b6c11cae71fdcb2be612885e9034f_83</lineID>
                phDoc(e) = tmpDoc;
helloWorld % <lineID version=octerine>46abfde9b1780a6feddb2a2416d85d46_84</lineID>
            end
helloWorld % <lineID version=octerine>e11262f65433794469efa9f38243b31a_85</lineID>
            JSON_string = savejson('cobDoc',phDoc);
helloWorld % <lineID version=octerine>df77127fc79013724dd7cdc81b075fd4_86</lineID>
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %% csv data file to disk
            %%%%%%%%%%%%%%%%%%%%%%%%%
            fileList = {};
helloWorld % <lineID version=octerine>0d9dcef21588ddde680707f02ce8ec08_87</lineID>
            fprintf(['starting save phase \n']);
helloWorld % <lineID version=octerine>9e38992f32149495e6ffe1a9d2084641_88</lineID>
            % save image
            fileList{end+1} = [oPath nm '_thumb.tif' ];
helloWorld % <lineID version=octerine>162111b7c31e95c0879849f1bce370f8_89</lineID>
            imwrite(thumb,fileList{end});
helloWorld % <lineID version=octerine>46dcc13caad789629717e001fdaca489_90</lineID>
            % save image
            fileList{end+1} = [oPath nm '_result.tif' ];
helloWorld % <lineID version=octerine>57fc917fd81ee9db0bce47f21e5c753d_91</lineID>
            saveas(h,fileList{end});
helloWorld % <lineID version=octerine>4354cbe8ac9c9389a4a9334b67f49522_92</lineID>
            % save mat file
            fileList{end+1} = [oPath nm '.mat'];
helloWorld % <lineID version=octerine>ef7b0d65865fea45e9efd53064ed6bde_93</lineID>
            save(fileList{end},'BB','fileName','S');
helloWorld % <lineID version=octerine>02bbecccf164d3f72dbd7644cb3fec52_94</lineID>
            % save the global parameters in file
            fileList{end+1} = [oPath nm '.csv'];
helloWorld % <lineID version=octerine>5dc03fc9db08e3d3fe5dc7acb5f59222_95</lineID>
            csvwrite(fileList{end},DATA);
helloWorld % <lineID version=octerine>58ea52bde022e3cb3b94bdbdbec6e7de_96</lineID>
            % save the width parameters
            fileList{end+1} = [oPath nm '_width_results.csv'];
helloWorld % <lineID version=octerine>044001c57f01ac926e53d46b10bb0482_97</lineID>
            csvwrite(fileList{end},S.widthProfile);
helloWorld % <lineID version=octerine>3f4e97a3e48b2802277e8efd17087199_98</lineID>
            % save the color values
            fileList{end+1} = [oPath nm '_cobRGB.csv'];
helloWorld % <lineID version=octerine>4e23c545387f41a41c90a07cda6b95c0_99</lineID>
            csvwrite(fileList{end},S.RGB);
helloWorld % <lineID version=octerine>97f88b54fbfd419b9afb4e5fb41c2036_100</lineID>
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % save JSON string
            %%%%%%%%%%%%%%%%%%%%%%%%%
            fileList{end+1} = [oPath nm '_jdoc.json'];
helloWorld % <lineID version=octerine>bccea8f12bc070eb5ba9b8125cb02a9c_101</lineID>
            fileID = fopen(fileList{end},'w');
helloWorld % <lineID version=octerine>9282ae052e3a1fa7bc1dda55987e5b3c_102</lineID>
            fprintf(fileID,strrep(JSON_string,'\/','\\/'));
helloWorld % <lineID version=octerine>54bdad7fda8751a0e14def6f0873af74_103</lineID>
            fclose(fileID);
helloWorld % <lineID version=octerine>0c9e2647b81280a5a7df558cbf244d4c_104</lineID>
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % push to iRODS
            %%%%%%%%%%%%%%%%%%%%%%%%%
            pushToiRods(rPath,fileList);
helloWorld % <lineID version=octerine>fbd21780845a4c16fe14552e4f3d8829_105</lineID>
            fprintf(['ending save phase \n']);
helloWorld % <lineID version=octerine>e5b7111e3429cdf77b32141d81517574_106</lineID>
        end
helloWorld % <lineID version=octerine>b07df4cea5b4d8c43da227186f6e0cd6_107</lineID>
        close all;
helloWorld % <lineID version=octerine>6e4c4aad7dac1d07a52d11cdf135920b_108</lineID>
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SAVE - end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %uncLog(['ph:l:'],{rawName},'set',algorithmName,algorithmVerion,{'5'},{'1'},1);
    catch ME
helloWorld % <lineID version=octerine>80b704acf1a62039e5c343167c7f4574_109</lineID>
        close all;
helloWorld % <lineID version=octerine>3a60424d5c9123e1672d681476dadc04_110</lineID>
        getReport(ME)
helloWorld % <lineID version=octerine>7eb5af479845e35543de8bf2899c8ef4_111</lineID>
        fprintf(['******error in:singleCobImage.m******\n']);
helloWorld % <lineID version=octerine>03f8f01a9fa109fda7aa22809bf86e6f_112</lineID>
        %uncLog(['ph:l:'],{rawName},'set',algorithmName,algorithmVerion,{'5'},{'2'},1);
    end
helloWorld % <lineID version=octerine>68e3e190b6f294ed0d96c016f0bc314f_113</lineID>
    close all
helloWorld % <lineID version=octerine>eba961e01cdf8554fedf6a7f9915e12f_114</lineID>
    fprintf(['****************************************************************************************\n']);
helloWorld % <lineID version=octerine>cbc6ef8e612e955d955254af07cdb254_115</lineID>
    fprintf(['Total Running Time: ' num2str(etime(clock,totalTimeInit)) '\n']);
helloWorld % <lineID version=octerine>755b4ce06554715254ef8d73fca96803_116</lineID>
    endString = 'Ending cob analysis algorithm. \n';
helloWorld % <lineID version=octerine>15e5a925fccdd19c425aead6c3625ee0_117</lineID>
    fprintf([endString,versionString]);
helloWorld % <lineID version=octerine>f674607bcfd575933223e40ea1dae480_118</lineID>
    fprintf(['****************************************************************************************\n']);
helloWorld % <lineID version=octerine>7101bdbb0c86e354499aee75972113e3_119</lineID>
    %uncLog(['ph:l:'],{rawName},'set',algorithmName,algorithmVerion,{'6'},{'1'},1);
end
helloWorld % <lineID version=octerine>829a5135b69188a0b7717694759cd6c8_120</lineID>


%{
<comBlock version=octerine>

</comBlock>
%}
