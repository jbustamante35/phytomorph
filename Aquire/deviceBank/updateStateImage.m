function [I] = updateStateImage(I,r,c,blockSZ,rgb,sat,subImage,AX)



    str1 = (r-1)*blockSZ(1)+1;
    str2 = (c-1)*blockSZ(2)+1;
    stp1 = str1 + blockSZ(1) - 1;
    stp2 = str2 + blockSZ(2) - 1;
    
    
    
    hsv = rgb2hsv(rgb);
    hsv(2) = sat;
    rgb = hsv2rgb(hsv);
    
    
    imageOut = makeColorBlock(blockSZ,rgb);
    
    if ~isempty(subImage)
        halphablend = vision.AlphaBlender;
        imageOut = step(halphablend,imageOut,subImage);
    end
    
    I(str1:stp1,str2:stp2,:) = imageOut;
    
    if ~isempty(AX)
        axes(AX);
        imshow(I);
    end
    
end