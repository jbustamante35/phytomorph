function [MASK,SKELETON,returnI,boundingBoxes] = cropPackAnalysisAndSave(cropPackage,GMModel,GMModel2)


    returnI = cropPackage.croppedImage;
    boundingBoxes = cropPackage.SEEDLING_BOX;

    for p = 1:numel(cropPackage.croppedImage)
        img = cropPackage.croppedImage{p};
        LAB = rgb2lab(img);
        sz = size(LAB);
        sam = LAB(:,:,1:3);
        sam = reshape(sam,[prod(sz(1:2)) 3]);



       
        kidx = GMModel.cluster(sam);
        kidx = reshape(kidx,sz(1:2));
        plantMask = kidx == 1;
        plantMask = bwareaopen(plantMask,500);
        fidx = find(plantMask);
        %out = flattenMaskOverlay(d{e}.croppedImage{p}/255,plantMask);


        kidx2 = GMModel2.cluster(sam(fidx,1:2));
        kidx(fidx) = kidx2+3;
        out = kidx==4;

        plantMask = out;
        plantMask = bwareaopen(plantMask,400);
        plantMask = imopen(plantMask,strel('disk',3,0));
        plantMask = imclearborder(plantMask);
        plantMask = bwareaopen(plantMask,400);
        [plantMask] = connectPlant2(plantMask);
        MASK{p} = plantMask;
        SKELETON{p} = bwmorph(plantMask,'skeleton',inf);
    end

end