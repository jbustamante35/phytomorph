function [] = gatherParametersForReader(imageStack,gatherType,dispMode)


    switch gatherType
        case 'histogramNormalize'
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% gather the normalization data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            toH = [];
            for e = 1:numel(imageStack)
                tmpH = imread(imageStack{e});
                toH = [toH;tmpH(:)];
                imshow(tmpH,[]);
                drawnow
                e
            end
    end


end