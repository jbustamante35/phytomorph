function [dataOut] = gatherParametersForReader(imageStack,gatherType,dispMode)
    toRand = true;
    dataOut = [];
    switch gatherType
        case 'histogram'
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% gather the normalization data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dataOut = [];
            % use only min of 100 or the number provided
            numToUse = min(100,numel(imageStack));
            if toRand;imageStack=imageStack(randperm(numel(imageStack)));end
            parfor e = 1:numToUse
                tmpH = imread(imageStack{e});
                dataOut = [dataOut;tmpH(:)];
                if dispMode
                    imshow(tmpH,[]);
                    drawnow
                end
                fprintf(['Done gathering histogram for:' num2str(e) '\n'])
            end
    end
end