function [S] = singleCobImage(fileName,noe,oPath,rPath,rawImage_scaleFactor,checkBlue_scaleFactor,defaultAreaPix,rho,addcut,baselineBlue,colRange1,colRange2,fill,toSave,toDisplay) % <lineID version=octerine>true</lineID>
varLogger(whos,1,1) % <lineID version=octerine>c5ba93e7751495989d54ad92e189c390_1</lineID>
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
varLogger(whos,2,1) % <lineID version=octerine>9893fad62374c668b9d8cbf41c7e9a90_2</lineID>
    algorithmVerion = '1.0';
varLogger(whos,3,1) % <lineID version=octerine>0469756e04a883714eac9bcf5bdc8b29_3</lineID>
    fprintf(['****************************************************************************************\n']);
varLogger(whos,4,1) % <lineID version=octerine>2a66c6ae90a13126afad49c020e9a6d5_4</lineID>
    versionString = 'Publication Version 1.0 - Monday, March 28, 2016. \n';
varLogger(whos,5,1) % <lineID version=octerine>dce6905ba216f85c841b37b1ebdad088_5</lineID>
    startString = 'Starting cob analysis algorithm. \n';
varLogger(whos,6,1) % <lineID version=octerine>8b2a24a703ea09030e385be7746f5ef7_6</lineID>
    fprintf([startString,versionString]);
varLogger(whos,7,1) % <lineID version=octerine>632e645bb04fe9323df0ba83553d3ddb_7</lineID>
    fprintf(['****************************************************************************************\n']);
varLogger(whos,8,1) % <lineID version=octerine>077a021de76929a9487d50cadbd0fae2_8</lineID>
    totalTimeInit = clock;
varLogger(whos,9,1) % <lineID version=octerine>aa10fbe6ece75f546a2158837852bcc9_9</lineID>
    try 
varLogger(whos,10,1) % <lineID version=octerine>8b99809b6c160b11e5de8bedb18aea4c_10</lineID>
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % INIT VARS - start
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['starting with variable and environment initialization.\n']);
varLogger(whos,11,1) % <lineID version=octerine>b40f16dc78d8177bcd427950d4dccc63_11</lineID>
        %%%%%%%%%%%%%%%%%%%%%%%
        % init reporting for %uncLog
        %%%%%%%%%%%%%%%%%%%%%%%
        [rawName ticket] = stripiTicket(fileName);
varLogger(whos,12,1) % <lineID version=octerine>14f80001ce625ec9e5330789f5a784f7_12</lineID>
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
varLogger(whos,13,1) % <lineID version=octerine>e24e1290f33dd8c5bf81d3d5694c332a_13</lineID>
        checkBlue_scaleFactor = StoN(checkBlue_scaleFactor);
varLogger(whos,14,1) % <lineID version=octerine>cc2a11c62ea84bef96402a3d3b0fe1d0_14</lineID>
        rawImage_scaleFactor = StoN(rawImage_scaleFactor);
varLogger(whos,15,1) % <lineID version=octerine>910417ce34917dd08180e9a29787117e_15</lineID>
        defaultAreaPix = StoN(defaultAreaPix);       
varLogger(whos,16,1) % <lineID version=octerine>26c3a1c6d5da619f90aee3708ab77c58_16</lineID>
        rho = StoN(rho);
varLogger(whos,17,1) % <lineID version=octerine>786db229590773d9c1f833658e73b12b_17</lineID>
        addcut = StoN(addcut);    
varLogger(whos,18,1) % <lineID version=octerine>c26bf630163ea2faae83bf20b4f3da7f_18</lineID>
        baselineBlue = StoN(baselineBlue);
varLogger(whos,19,1) % <lineID version=octerine>bd2def8301218f701c5a204683684de7_19</lineID>
        colRange1 = StoN(colRange1);
varLogger(whos,20,1) % <lineID version=octerine>fb1f6a9e3e45b256232fa56156697791_20</lineID>
        colRange2 = StoN(colRange2);
varLogger(whos,21,1) % <lineID version=octerine>592545dbeb9530640e8b60be8c8bf344_21</lineID>
        fill = StoN(fill);
varLogger(whos,22,1) % <lineID version=octerine>2008458e7521638ccb6242074588d4e1_22</lineID>
        %%%%%%%%%%%%%%%%%%%%%%%
        % print out the input variables
        %%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['FileName:' fileName '\n']);
varLogger(whos,23,1) % <lineID version=octerine>172cccaa1fd4a4d517cbdaf203ae1cbd_23</lineID>
        fprintf(['Number of Ears:' num2str(noe) '\n']);
varLogger(whos,24,1) % <lineID version=octerine>6197a99bc6432c8e9e414898e4cbdf49_24</lineID>
        fprintf(['OutPath:' oPath '\n']);     
varLogger(whos,25,1) % <lineID version=octerine>598d63745b18d81ad90f731803ef384a_25</lineID>
        fprintf(['Image resize in checkBlue:' num2str(checkBlue_scaleFactor) '\n']); 
varLogger(whos,26,1) % <lineID version=octerine>0b7ec9f0ba1ecbef4fdcaca0b104484f_26</lineID>
        fprintf(['Raw image resize:' num2str(rawImage_scaleFactor) '\n']);  
varLogger(whos,27,1) % <lineID version=octerine>43c0f529bd722d85c90795f42a680d78_27</lineID>
        fprintf(['Threshold noise size:' num2str(defaultAreaPix) '\n']);
varLogger(whos,28,1) % <lineID version=octerine>9393dd324369f3ee2ec95e6df4d68c98_28</lineID>
        fprintf(['The radius of color circle:' num2str(rho) '\n']);
varLogger(whos,29,1) % <lineID version=octerine>9820988ca100a6d7d97cbaf1943eaae7_29</lineID>
        fprintf(['The boarder handle for checkBlue:' num2str(addcut) '\n']);
varLogger(whos,30,1) % <lineID version=octerine>650e1805c5fcbf089eef0d01b83ee324_30</lineID>
        fprintf(['Baseline threshold to remove blue header:' num2str(baselineBlue) '\n']);
varLogger(whos,31,1) % <lineID version=octerine>c048d5bea5e854d547ceaf30ae974d57_31</lineID>
        fprintf(['Background Color Range I:' num2str(colRange1) '\n']);
varLogger(whos,32,1) % <lineID version=octerine>2d4d1fc20db6dc1efc356d7fc29e7ba1_32</lineID>
        fprintf(['Background Color Range II:' num2str(colRange2) '\n']);
varLogger(whos,33,1) % <lineID version=octerine>848baf9792d644506d48dde0591a3997_33</lineID>
        fprintf(['The radius of disk for closing:' num2str(fill) '\n']);
varLogger(whos,34,1) % <lineID version=octerine>8fbc9b693a1afec1103e3d27a9a82bc6_34</lineID>
        %%%%%%%%%%%%%%%%%%%%%%%
        % make output directory
        %%%%%%%%%%%%%%%%%%%%%%%
        mkdir(oPath);
varLogger(whos,35,1) % <lineID version=octerine>53563fc024df0a39670f599133b90c45_35</lineID>
        [pth nm ext] = fileparts(fileName);
varLogger(whos,36,1) % <lineID version=octerine>a711e45f09c459adc0e73fc54a782668_36</lineID>
        fprintf(['ending with variable and environment initialization.\n']);
varLogger(whos,37,1) % <lineID version=octerine>35f39ebb18b434e8f90720f61794617f_37</lineID>
        %%%%%%%%%%%%%%%%%%%%%%%
        % read the image and take off the blue strip for bar code
        %%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['starting with image load.\n']);
varLogger(whos,38,1) % <lineID version=octerine>d5bdd9bc0dd1ea67bd61ee4dc3fd1dfb_38</lineID>
        %uncLog(['ph:l:'],{rawName},'set',algorithmName,algorithmVerion,{'3'},{'1'},1);
        I = imread(fileName);
varLogger(whos,39,1) % <lineID version=octerine>34bb12205752424245a32122b8e881db_39</lineID>
        % remove 4th pane for some images
        I = I(:,:,1:3);
varLogger(whos,40,1) % <lineID version=octerine>9c06a020897da2bf238445b61acf351c_40</lineID>
        %uncLog(['ph:l:'],{rawName},'set',algorithmName,algorithmVerion,{'4'},{'1'},1);
        thumb = imresize(I,.25);
varLogger(whos,41,1) % <lineID version=octerine>8eea5ad80f6143d5de20373bd1063e80_41</lineID>
        % rawImage_scaleFactor to lower 'DPI' effect, by fraction
        % If resize factor is 1, do not excecute imresize
        if rawImage_scaleFactor ~= 1;I = imresize(I,rawImage_scaleFactor);end
varLogger(whos,42,1) % <lineID version=octerine>4877852744964fa14432cd9d515edf2d_42</lineID>
        % check blue header and remove
        I = checkBlue(I,checkBlue_scaleFactor,addcut,baselineBlue);
varLogger(whos,43,1) % <lineID version=octerine>41d56710676f1e05e77d0b53d8e86c3e_43</lineID>
        fprintf(['ending with image load.\n']);
varLogger(whos,44,1) % <lineID version=octerine>553bd2ff7bb2faac792875c3facecb3f_44</lineID>
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % INIT VARS - end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ANALYSIS - start
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
        % main measurement code call
        fprintf(['starting with image analysis \n']);
varLogger(whos,45,1) % <lineID version=octerine>ce074b48b91ca09c992a75386f5d9ca4_45</lineID>
        % convert to single
        %I = single(I)/255;
        % run the analysis
        [BB S] = maizeCob(I,noe,defaultAreaPix,colRange1,colRange2,fill);
varLogger(whos,46,1) % <lineID version=octerine>6fdbc48390cce5db3ac1b6866819d14a_46</lineID>
        % stack the results from the bounding box
        DATA = [];
varLogger(whos,47,1) % <lineID version=octerine>0896d9ebc3a5908a52f72fd60868f0d8_47</lineID>
        for b = 1:numel(BB)                
varLogger(whos,48,1) % <lineID version=octerine>325993baa053888d8e17911ecd84e9fb_48</lineID>
            DATA = [DATA;[BB{b}(3:4) S.average_WIDTH(b)]];        
varLogger(whos,49,1) % <lineID version=octerine>22d9dceee375e2a7fda7910151598b1b_49</lineID>
        end
varLogger(whos,50,1) % <lineID version=octerine>f542f5ee1a7f40bdc44f43c5070a67b5_50</lineID>
        uDATA = mean(DATA,1);
varLogger(whos,51,1) % <lineID version=octerine>44d86ddb0e6d2001a3f4e57647b1c2f8_51</lineID>
        sDATA = std(DATA,1,1);
varLogger(whos,52,1) % <lineID version=octerine>646cf271ab82f62c757104edebe82d9b_52</lineID>
        DATA = reshape(DATA',[1 numel(DATA)]);
varLogger(whos,53,1) % <lineID version=octerine>ecd41f00ca630db3c6ac3a44e66536d7_53</lineID>
        fprintf(['ending with image analysis \n']);
varLogger(whos,54,1) % <lineID version=octerine>46bd7a4cf3d49346b1cf8cb8eeac60da_54</lineID>
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ANALYSIS - start
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DISLAY - start
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
        if toDisplay
varLogger(whos,55,1) % <lineID version=octerine>9d0a2e2d2e86ed0ba445a563fbc15bd7_55</lineID>
            fprintf(['starting with image and results display \n']);
varLogger(whos,56,1) % <lineID version=octerine>fb05075f039e87558b735aa3c00ca547_56</lineID>
            h = image(I);
varLogger(whos,57,1) % <lineID version=octerine>970830f708c2914a71ba2212c9964daf_57</lineID>
            hold on
varLogger(whos,58,1) % <lineID version=octerine>b55aa0b87df5a80cea8145157cb7d06c_58</lineID>
            % make rectangle around the cob
            for b = 1:numel(BB)
varLogger(whos,59,1) % <lineID version=octerine>3db420eb01c01d031b6f37fde6ce0c16_59</lineID>
                rectangle('Position',BB{b},'EdgeColor','r');
varLogger(whos,60,1) % <lineID version=octerine>22f7924b6a3cb714345f7ae1280b9843_60</lineID>
                UR = BB{b}(1:2);
varLogger(whos,61,1) % <lineID version=octerine>86569aa24a3f476de7981feea4dcc428_61</lineID>
                UR(1) = UR(1) + BB{b}(3);
varLogger(whos,62,1) % <lineID version=octerine>b6d0983af133710a5873061dfb00070d_62</lineID>
                rectangle('Position',[UR rho rho],'EdgeColor','none','Curvature',[1 1],'FaceColor',S.RGB(b,:)/255);
varLogger(whos,63,1) % <lineID version=octerine>d6fb5742237b5fae995084094a905b38_63</lineID>
            end
varLogger(whos,64,1) % <lineID version=octerine>e550b9200777416a604e962012031f28_64</lineID>
            %%%%%%%%%%%%%%%%%%%%%%%
            % format the image
            axis equal;axis off;drawnow;set(gca,'Position',[0 0 1 1]);
varLogger(whos,65,1) % <lineID version=octerine>c66e53da31092ccde2aec49f6b327d17_65</lineID>
            fprintf(['ending with image and results display \n']);
varLogger(whos,66,1) % <lineID version=octerine>f5664b7bece9f4c94e07b48b395f2994_66</lineID>
        end
varLogger(whos,67,1) % <lineID version=octerine>153c965eaeaea391a658f3c122977f69_67</lineID>
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DISLAY - end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SAVE - start
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
        if toSave
varLogger(whos,68,1) % <lineID version=octerine>d03148f01dec3b80c68afd390c9c8e21_68</lineID>
            %% data to database
            for e = 1:numel(BB)
varLogger(whos,69,1) % <lineID version=octerine>137d6efc27a0d4cc1f226cb8b1ab4913_69</lineID>
                %uncLog(['ph:p:'],{rawName},'set',algorithmName,[algorithmVerion],{['cob:length:' num2str(e)]},{num2str(BB{b}(4))},1);
                %uncLog(['ph:p:'],{rawName},'set',algorithmName,[algorithmVerion],{['cob:width:' num2str(e)]},{num2str(BB{b}(3))},1);
            end
varLogger(whos,70,1) % <lineID version=octerine>3f4f467e87d295dd22d101a9bb3ce11b_70</lineID>
            %uncLog(['ph:p:'],{rawName},'set',algorithmName,[algorithmVerion],{['cob:average_length:' num2str(e)]},{num2str(uDATA(2))},1);
            %uncLog(['ph:p:'],{rawName},'set',algorithmName,[algorithmVerion],{['cob:std_length:' num2str(e)]},{num2str(sDATA(2))},1);
            %uncLog(['ph:p:'],{rawName},'set',algorithmName,[algorithmVerion],{['cob:std_width:' num2str(e)]},{num2str(sDATA(1))},1);
            %uncLog(['ph:p:'],{rawName},'set',algorithmName,[algorithmVerion],{['cob:std_width:' num2str(e)]},{num2str(sDATA(1))},1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %% spin-up JSON output format
            %%%%%%%%%%%%%%%%%%%%%%%%%
            linkTo = stripiTicket(fileName);
varLogger(whos,71,1) % <lineID version=octerine>73c7d2a3c81bf89182e3e60c1b32ee96_71</lineID>
            linkPath = stripiTicket(rPath);
varLogger(whos,72,1) % <lineID version=octerine>7ab7b4a55dfd78aee2566aadb96027cd_72</lineID>
            [jP,tN] = fileparts(fileName);
varLogger(whos,73,1) % <lineID version=octerine>ffadc09e9f73b4d3134b6b908ed72e2a_73</lineID>
            for e = 1:numel(BB)
varLogger(whos,74,1) % <lineID version=octerine>b1f4e6a319c4cb74c17811f754544800_74</lineID>
                tmpDoc = [];
varLogger(whos,75,1) % <lineID version=octerine>448db5b7f3c8a4cc0d3996e3cfc5e153_75</lineID>
                tmpDoc = generatePhenotypeNode(tmpDoc,linkTo,'orginalImage','orginalImage');
varLogger(whos,76,1) % <lineID version=octerine>c1733af4b539e952772657e2479fba80_76</lineID>
                tmpDoc = generatePhenotypeNode(tmpDoc,[linkPath filesep tN '_thumb.tif' ],'thumbNail','thumbNail');
varLogger(whos,77,1) % <lineID version=octerine>bf5a144310c1bc053fb48baf635789bf_77</lineID>
                tmpDoc = generatePhenotypeNode(tmpDoc,S.widthProfile(e,:),{'along','position'},'widthProfile');
varLogger(whos,78,1) % <lineID version=octerine>e1737f3badcdc517e7ad178aa6a5d4b3_78</lineID>
                tmpDoc = generatePhenotypeNode(tmpDoc,S.RGB(e,:),{'along','colorIndex'},'rgb_color');
varLogger(whos,79,1) % <lineID version=octerine>a5a96b0d31db0451f8fdf4e318c1b562_79</lineID>
                tmpDoc = generatePhenotypeNode(tmpDoc,BB{e}(3),{'maxWidth'},'maxWidth');
varLogger(whos,80,1) % <lineID version=octerine>4b4906ea679737d506f8422060128d3f_80</lineID>
                tmpDoc = generatePhenotypeNode(tmpDoc,S.average_WIDTH(e),{'averageWidth'},'averageWidth');
varLogger(whos,81,1) % <lineID version=octerine>854fd445b48f970faf02f343d1c4c511_81</lineID>
                tmpDoc = generatePhenotypeNode(tmpDoc,BB{e}(4),{'length'},'length');
varLogger(whos,82,1) % <lineID version=octerine>3ca412d91f2e712a4ca6a6d01f064c0d_82</lineID>
                tmpDoc = generatePhenotypeNode(tmpDoc,BB{e},{'boundingBox'},'boundingBox');
varLogger(whos,83,1) % <lineID version=octerine>f137fe33b965513140a319091e99f0c8_83</lineID>
                phDoc(e) = tmpDoc;
varLogger(whos,84,1) % <lineID version=octerine>45815b42dbb6bf236597194c4f6e3b45_84</lineID>
            end
varLogger(whos,85,1) % <lineID version=octerine>9ecc27671cf0f6bb1969f45793648c1b_85</lineID>
            JSON_string = savejson('cobDoc',phDoc);
varLogger(whos,86,1) % <lineID version=octerine>cbf722a6969e5ef00db7e105ebab7220_86</lineID>
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %% csv data file to disk
            %%%%%%%%%%%%%%%%%%%%%%%%%
            fileList = {};
varLogger(whos,87,1) % <lineID version=octerine>7d886d3360b96db5505997ebfbdcbdd7_87</lineID>
            fprintf(['starting save phase \n']);
varLogger(whos,88,1) % <lineID version=octerine>05ff79999bfe8fe0186ae2abf2e78601_88</lineID>
            % save image
            fileList{end+1} = [oPath nm '_thumb.tif' ];
varLogger(whos,89,1) % <lineID version=octerine>0b826025a7cdfd9ca4dbe851095947f6_89</lineID>
            imwrite(thumb,fileList{end});
varLogger(whos,90,1) % <lineID version=octerine>614fb4d011172127d0ade53b97c53dbc_90</lineID>
            % save image
            fileList{end+1} = [oPath nm '_result.tif' ];
varLogger(whos,91,1) % <lineID version=octerine>4df3ac2613b8064e6311ef6b6cacb316_91</lineID>
            saveas(h,fileList{end});
varLogger(whos,92,1) % <lineID version=octerine>98a72e3440d737028b22d7fbde206a50_92</lineID>
            % save mat file
            fileList{end+1} = [oPath nm '.mat'];
varLogger(whos,93,1) % <lineID version=octerine>19b3458102c497a899b3a3c90ccbb619_93</lineID>
            save(fileList{end},'BB','fileName','S');
varLogger(whos,94,1) % <lineID version=octerine>944a2a5633195ec711cce1d9df83cbf1_94</lineID>
            % save the global parameters in file
            fileList{end+1} = [oPath nm '.csv'];
varLogger(whos,95,1) % <lineID version=octerine>31f157a278d3940ccdce5380dc1b14fa_95</lineID>
            csvwrite(fileList{end},DATA);
varLogger(whos,96,1) % <lineID version=octerine>8d3af5e58df6a47da2f5a0c377ef4df4_96</lineID>
            % save the width parameters
            fileList{end+1} = [oPath nm '_width_results.csv'];
varLogger(whos,97,1) % <lineID version=octerine>25ed58468057442480c259439d4d46ce_97</lineID>
            csvwrite(fileList{end},S.widthProfile);
varLogger(whos,98,1) % <lineID version=octerine>f37aee1ffb0d53aa5d3a61fc4f1e136f_98</lineID>
            % save the color values
            fileList{end+1} = [oPath nm '_cobRGB.csv'];
varLogger(whos,99,1) % <lineID version=octerine>5652ef590a103f3be5380a602829d4f0_99</lineID>
            csvwrite(fileList{end},S.RGB);
varLogger(whos,100,1) % <lineID version=octerine>d3f761eda316e0e0fedf58a147b1e7f6_100</lineID>
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % save JSON string
            %%%%%%%%%%%%%%%%%%%%%%%%%
            fileList{end+1} = [oPath nm '_jdoc.json'];
varLogger(whos,101,1) % <lineID version=octerine>07d47a3fc06c09842f6627eae4cd72a1_101</lineID>
            fileID = fopen(fileList{end},'w');
varLogger(whos,102,1) % <lineID version=octerine>092b6223e7148969464814aa72e148ef_102</lineID>
            fprintf(fileID,strrep(JSON_string,'\/','\\/'));
varLogger(whos,103,1) % <lineID version=octerine>4f3e6454ede68b2fab4e145f396bad12_103</lineID>
            fclose(fileID);
varLogger(whos,104,1) % <lineID version=octerine>51da500075c2ed87e6f5643b6b8828eb_104</lineID>
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % push to iRODS
            %%%%%%%%%%%%%%%%%%%%%%%%%
            pushToiRods(rPath,fileList);
varLogger(whos,105,1) % <lineID version=octerine>caa82e6311293e7810a937069b1721ee_105</lineID>
            fprintf(['ending save phase \n']);
varLogger(whos,106,1) % <lineID version=octerine>7f07252d1e245d7f54203244e68cc96b_106</lineID>
        end
varLogger(whos,107,1) % <lineID version=octerine>0dd68832e5ac02d97920183d234ece77_107</lineID>
        close all;
varLogger(whos,108,1) % <lineID version=octerine>70b45526223341ca80fb21f793f144cf_108</lineID>
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SAVE - end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %uncLog(['ph:l:'],{rawName},'set',algorithmName,algorithmVerion,{'5'},{'1'},1);
    catch ME
varLogger(whos,109,1) % <lineID version=octerine>317e9be9705a3ba6bfc62aa1052abd19_109</lineID>
        close all;
varLogger(whos,110,1) % <lineID version=octerine>db634ce0c69346ae5ca3dbf9af27cf91_110</lineID>
        getReport(ME)
varLogger(whos,111,1) % <lineID version=octerine>daea53e26da734e4dcf96cfecd797acc_111</lineID>
        fprintf(['******error in:singleCobImage.m******\n']);
varLogger(whos,112,1) % <lineID version=octerine>1f88065dceb973dd95a4bf8219cc5655_112</lineID>
        %uncLog(['ph:l:'],{rawName},'set',algorithmName,algorithmVerion,{'5'},{'2'},1);
    end
varLogger(whos,113,1) % <lineID version=octerine>f97727d99005a17a3bda396c10ebb6ae_113</lineID>
    close all
varLogger(whos,114,1) % <lineID version=octerine>5c730e2d8e5da3cac657eecf65894fe6_114</lineID>
    fprintf(['****************************************************************************************\n']);
varLogger(whos,115,1) % <lineID version=octerine>ce89a464a5e6ca39285701fdf04065fb_115</lineID>
    fprintf(['Total Running Time: ' num2str(etime(clock,totalTimeInit)) '\n']);
varLogger(whos,116,1) % <lineID version=octerine>2530f48059e0e9f62d14897d06e861cd_116</lineID>
    endString = 'Ending cob analysis algorithm. \n';
varLogger(whos,117,1) % <lineID version=octerine>fbce7a1912152870a59a4000d6748b69_117</lineID>
    fprintf([endString,versionString]);
varLogger(whos,118,1) % <lineID version=octerine>a9db003d8adda087885c13d58789816a_118</lineID>
    fprintf(['****************************************************************************************\n']);
varLogger(whos,119,1) % <lineID version=octerine>edd8526a0d27303969d04394879483d7_119</lineID>
    %uncLog(['ph:l:'],{rawName},'set',algorithmName,algorithmVerion,{'6'},{'1'},1);
end
varLogger(whos,120,1) % <lineID version=octerine>0aa5cf37277fe5da5db482f80441bd5f_120</lineID>

