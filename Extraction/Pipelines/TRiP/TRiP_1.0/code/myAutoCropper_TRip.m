function [] = myAutoCropper_TRip(FileList,disp,oPath)

    %%%%%%%%%%%%%%%%%%%%%%%%%
    oPath = [oPath];
    I = [];
    for e = 1:1:numel(FileList)
        tmp = imread(FileList{e});
        Lab = double(rgb2lab(tmp))/255;
        tmp = imfilter(Lab(:,:,2),fspecial('gaussian',[31 31],11),'replicate');
        I(:,:,e) = tmp;
        %imshow(I(:,:,e),[]);
        fprintf(['Done loading :' num2str(e) ':' num2str(numel(FileList)) '\n'])
    end
    MSK = bindVec(-std(I,1,3).*mean(I,3));
    MSK = MSK > graythresh(MSK);
    MSK = bwlarge(MSK,9);
    R = regionprops(MSK,'Centroid');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    for e = 1:1:numel(FileList)
        tmp = imread(FileList{e});
        [pth,nm,ext] = fileparts(FileList{e});
        if disp
            imshow(tmp,[]);
            hold on
            for r = 1:numel(R)
                box = R(r).Centroid - 50;
                rectangle('Position',[box 100 100],'EdgeColor','r')
            end
            drawnow
            hold off
        end
        
        for r = 1:numel(R)
            tmp_oPath = [oPath 'crop_plant' num2str(r) filesep];
            mkdir(tmp_oPath)
            box = R(r).Centroid - 50;
            cropped = imcrop(tmp,[box 100 100]);
            imwrite(cropped,[tmp_oPath 'crop_' nm '.tif']);
        end
        fprintf(['Done crop :' num2str(e) ':' num2str(numel(FileList)) '\n'])
    end
    
    
end

%{
    
   
%}

