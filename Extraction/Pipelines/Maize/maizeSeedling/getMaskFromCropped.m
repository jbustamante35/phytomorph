function [MASK,SKELETON,returnI,boundingBoxes] = getMaskFromCropped(cropPackage,E,U,GMModel_Color)


    returnI = cropPackage.croppedImage;
    boundingBoxes = cropPackage.SEEDLING_BOX;

    for e = 1:numel(returnI)

        I = rgb2hsv(returnI{e});
        Io = returnI{e};
        I = I(:,:,3);
        sig = mean(I,2);
        sig = sig((end-200):end);
        sig = gradient(sig);
        sig = im2col(sig,[7 1],'sliding');
        tC = PCA_REPROJ_T(sig,E,U);
        tC = tC(1:180);
        [~,midx] = min(tC);
        midx = size(I,1) - (200 - midx);
        Ic = Io(1:midx,:,:);
        returnI{e} = Ic;
        Lab = rgb2lab(Ic);
        sz = size(Lab);
        CL = reshape(Lab(:,:,2:3),[prod(sz(1:2)) 2]);



        kidx = GMModel_Color.cluster(CL);
        kidx = reshape(kidx,sz(1:2));
        plantMask = kidx == 2;
        plantMask = bwareaopen(plantMask,300);
        
        [plantMask] = connectPlant2(plantMask);
        MASK{e} = plantMask;
        SKELETON{e} = bwmorph(plantMask,'skeleton',inf);
        
        
    end
end