function [] = singulatedKernels(fileName,oPath)
    try
        mkdir(oPath);
        [pth,nm,ext] = fileparts(fileName);
        % read image
        I = double(imread(fileName));
        % make gray scale color image
        if size(I,3) == 1
            I = cat(3,I,I,I);
        end
        % get centers
        [C,imageSZ,RxC,imageSZ,BOX] = getCenters_singulated(fileName);
        % align centers
        C = realignCenters_singulated(fileName,C,BOX);
        % get row count
        rowCount = RxC(1);


        h1 = figure;
        image(I/255);
        axis equal
        hold on
        plot(C(:,2),C(:,1),'g*')

        csvOut = {'Area','Width','Height','Eccentricity','Perimeter'};


        for e = 1:size(C,1)
            tmpBOX = [fliplr(C(e,:))-BOX(3:4) 2*BOX(3:4)];
            rectangle('Position',tmpBOX,'EdgeColor','r');
            tmpKernel = imcrop(I,tmpBOX);



            [level,tmpMASK] = getThresholdLevel_sigulated(tmpKernel);
            
            % threshold the image
            tmpMASK = single(tmpMASK) > level;
            tmpMASK = imclearborder(tmpMASK);
            tmpMASK = bwlarge(tmpMASK);

            
            tmpKernel = flattenMaskOverlay(tmpKernel/255,tmpMASK,.6,'r');
            if any(tmpMASK(:)==1)
                H(e) = sum(any(tmpMASK,2));
                W(e) = sum(any(tmpMASK,1));
                A(e) = sum(tmpMASK(:));
                
                
                R = regionprops(logical(tmpMASK),'Eccentricity','Perimeter','Centroid');

                
            else
                H(e) = NaN;
                W(e) = NaN;
                A(e) = NaN;
                R(1).Eccentricity = NaN;
                R(1).Perimeter = NaN;
                cp = size(tmpKernel)/2;
                cp = fliplr(cp(1:2));
                R(1).Centroid =  cp;
            end


            csvOut{e+1,1} = A(e);
            csvOut{e+1,2} = W(e);
            csvOut{e+1,3} = H(e);
            csvOut{e+1,4} = R(1).Eccentricity;
            csvOut{e+1,5} = R(1).Perimeter;



            tmpV = sum(tmpMASK,2);
            shapeMatrix(e,:) = interp1(1:numel(tmpV),tmpV,linspace(1,numel(tmpV),1000));


            figure(h1);
            textPOS = fliplr(C(e,:))-BOX(3:4) - [20 20];
            text(textPOS(1),textPOS(2),num2str(e),'BackgroundColor','w')



            h2 = figure;
            image(tmpKernel);
            axis equal
            hold on

            cp = size(tmpKernel)/2;
            cp = fliplr(cp(1:2));
            cp = R(1).Centroid;
            VERT = [[cp(1) cp(2)-H(e)/2];[cp(1) cp(2)+H(e)/2]];
            plot(VERT(:,1),VERT(:,2),'g')

            HOR = [[cp(1)-W(e)/2 cp(2)];[cp(1)+W(e)/2 cp(2)]];
            plot(HOR(:,1),HOR(:,2),'g')
            hold off
            saveas(h2,[oPath nm '_' num2str(e) '.jpg']);
            imwrite(tmpKernel,[oPath nm '_' num2str(e) '_raw.tif']);
            
            %waitforbuttonpress
            close(h2);

            



            drawnow
        end

        saveas(h1,[oPath nm '_raster.jpg']);
        cell2csv([oPath nm '_kernelMeasurements.csv'],csvOut);
        csvwrite([oPath nm '_widthProfiles.csv'],shapeMatrix);

        close all
    catch ME
        getReport(ME)
    end
end

%{
    fileName = '~/Download/1.tiff';
   
    fileName = '~/test.tif';
    fileName = '~/Download/Plates1-2_abgerminal.tiff';
    fileName = '~/Download/Plates1-2_abgerminalCROPPED.tiff';
    %fileName = '~/Download/Plates11-12_longitudinalCROPPED.tiff';
    %fileName = '~/Download/Plates3-4_capsCROPPED.tiff';
    fileName = '~/Plates5-6_abgerminalCROPPED.tiff';
    fileName = '~/Plates7-8_longitudinalCROPPED.tiff';
    fileName = '~/Mask of Plates1-2_germinal_rotated_glumeless - Copy.tiff';
    oPath = './output/';
    singulatedKernels(fileName,oPath);

    

%}