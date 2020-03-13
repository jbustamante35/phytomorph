function [Ir,other] = steveAlign(I,other)
    try
        ptsOriginal  = detectSURFFeatures(I(:,:,1));
        ptsDistorted = detectSURFFeatures(I(:,:,2));


        [featuresOriginal,validPtsOriginal] = ...
            extractFeatures(I(:,:,1),ptsOriginal);
        [featuresDistorted,validPtsDistorted] = ...
            extractFeatures(I(:,:,2),ptsDistorted);
        index_pairs = matchFeatures(featuresOriginal,featuresDistorted);
        matchedPtsOriginal  = validPtsOriginal(index_pairs(:,1));
        matchedPtsDistorted = validPtsDistorted(index_pairs(:,2));
        %{
        figure; 
        showMatchedFeatures(I(:,:,1),I(:,:,2),...
            matchedPtsOriginal,matchedPtsDistorted);
        %}
        [tform,inlierPtsDistorted,inlierPtsOriginal] = ...
            estimateGeometricTransform(matchedPtsDistorted,matchedPtsOriginal,'affine');

        %{
        figure; 
        showMatchedFeatures(I(:,:,1),I(:,:,2),...
            inlierPtsOriginal,inlierPtsDistorted);
        title('Matched inlier points');
        %}

        outputView = imref2d(size(I(:,:,1)));
        Ir = imwarp(I(:,:,2),tform,'OutputView',outputView);
        
        other = imwarp(other,tform,'OutputView',outputView);
        %{
        figure; imshow(Ir); 
        title('Recovered image');
        %}
    catch ME
        ME
    end
end
    