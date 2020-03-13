function [percentGreen,uNDVI,sNDVI] = singleGreenFieldImage(fileName,GMM,GMM_NIR,type,oPath)
    warning off
    mkdir(oPath);
    fileList = {};
    [pth,nm,ext] = fileparts(fileName);

    
    percentGreen = 0;
    uNDVI = 0;
    sNDVI = 0;
    if type == 'RGB'

        I = imread(fileName);
        % convert to Lab
        Lab = rgb2lab(I);
        
        I = double(I);


        if exist('GMM')
            sz = size(I);
            tmp = double(reshape(Lab,[prod(sz(1:2)) sz(3)]));
            kidx = GMM.cluster(tmp);
            kidx = reshape(kidx,sz(1:2));
            plantMask = kidx == 2;
            out = flattenMaskOverlay(double(I)/255,plantMask);

            percentGreen = sum(plantMask(:)) / numel(plantMask);

            DISKMask = Lab(:,:,1) > 99;
            DISKMask = bwlarge(DISKMask);
            out = flattenMaskOverlay(out,DISKMask,.5,'b');
            if ~isdeployed()
                %imshow(out,[]);
            end
            imwrite(out,[oPath nm '_returnOverlayRGB' '.jpg']);

            for k = 1:size(I,3)
                tmp = I(:,:,k);
                HISTO(:,k) = hist(tmp(find(plantMask))/255,linspace(0,1,256));
            end
        end


        csvwrite([oPath nm '_return_RGB_CSV' '.csv'],[percentGreen]);
        
        


        tmpDoc = [];
        tmpDoc = generatePhenotypeNode(tmpDoc,nm,'FileName','FileName');
        tmpDoc = generatePhenotypeNode(tmpDoc,percentGreen,'percentGreen','percentGreen');
        tmpDoc = generatePhenotypeNode(tmpDoc,HISTO,{'RGBcolorValue','channel'},'RGB_histogramm');
        
        JSON_string = savejson('rgbDoc',tmpDoc);
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % save JSON string
        %%%%%%%%%%%%%%%%%%%%%%%%%
        fileList{end+1} = [oPath nm '_jdoc.json'];
        fileID = fopen(fileList{end},'w');
        fprintf(fileID,strrep(JSON_string,'\/','\\/'));
        fclose(fileID);
        %%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
    else
        [pth,nm,ext] = fileparts(fileName);
        I = imread(fileName);
        Lab = rgb2lab(I);

        I = double(I);

        if exist('GMM_NIR')
            sz = size(I);
            tmp = double(reshape(Lab,[prod(sz(1:2)) sz(3)]));
            kidx = GMM_NIR.cluster(tmp);
            kidx = reshape(kidx,sz(1:2));
            plantMask = kidx == 2;
            out = flattenMaskOverlay(double(I)/255,plantMask);



            DISKMask = Lab(:,:,1) > 98;
            DISKMask = bwlarge(DISKMask);
            DISKMask = imerode(DISKMask,strel('disk',2,0));

            for k = 1:3
                tmp = I(:,:,k);
                urgb(k) = mean(tmp(find(DISKMask)));
            end
            delta1 = 255 - urgb(1);
            I(:,:,1) = I(:,:,1) + delta1;


            NDVI = (I(:,:,1) - I(:,:,2)).*(I(:,:,1) + I(:,:,2)).^-1;

            uNDVI = mean(NDVI(find(plantMask)));
            sNDVI = std(NDVI(find(plantMask)));

            out = flattenMaskOverlay(out,DISKMask,.5,'b');
            
            if ~isdeployed()
                %imshow(out,[]);
            end
            
            
            imwrite(out,[oPath nm '_returnOverlayNIR' '.jpg']);
            
            csvwrite([oPath nm '_return_NIR_CSV' '.csv'],[uNDVI sNDVI]);
            
            
            
            
            HISTO = hist(NDVI(find(plantMask)),linspace(0,1,256));
            
            
            tmpDoc = [];
            tmpDoc = generatePhenotypeNode(tmpDoc,nm,'FileName','FileName');
            tmpDoc = generatePhenotypeNode(tmpDoc,uNDVI,'AverageNDVI','AverageNDVI');
            tmpDoc = generatePhenotypeNode(tmpDoc,sNDVI,'StdNDVI','StdNDVI');
            tmpDoc = generatePhenotypeNode(tmpDoc,HISTO',{'NDVIcolorValue'},'NDVI_histogramm');
            JSON_string = savejson('ndviDoc',tmpDoc);
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % save JSON string
            %%%%%%%%%%%%%%%%%%%%%%%%%
            fileList{end+1} = [oPath nm '_jdoc.json'];
            fileID = fopen(fileList{end},'w');
            fprintf(fileID,strrep(JSON_string,'\/','\\/'));
            fclose(fileID);
            %%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        
    end
end