function [valError,cons,trainedNet] = makeObjFcn_classify(...
                                    optVars,...
                                    trainImages,...
                                    trainLabels,...
                                    valImages,...
                                    valLabels,...
                                    toSLOW,...
                                    DISP,...
                                    MaxE,...
                                    exeEnvironment,...
                                    nDIMS,...
                                    TYPE)

        if nargin == 10
            TYPE = 'sequence';
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % create options for major train
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        options = trainingOptions('sgdm',...
            'InitialLearnRate',optVars.InitialLearnRate,...
            'Momentum',optVars.Momentum,...
            'MaxEpochs',MaxE,...
            'L2Regularization',optVars.L2Regularization,...
            'Shuffle','every-epoch',...
            'Verbose',true,...
            'ExecutionEnvironment',exeEnvironment,...
            'Plots',DISP);




        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % make layers
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        layers = [ ...
            sequenceInputLayer(nDIMS)
            lstmLayer(optVars.lstmStates,'OutputMode',TYPE)
            fullyConnectedLayer(2)
            softmaxLayer
            classificationLayer];



        trainedNet = trainNetwork(trainImages,trainLabels,layers,options);



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % create options for slow train
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if slow then get steady state ~ for stuff
        if toSLOW
            options = trainingOptions('sgdm',...
            'InitialLearnRate',optVars.InitialLearnRate/10,...
            'Momentum',optVars.Momentum,...
            'MaxEpochs',MaxE,...
            'L2Regularization',optVars.L2Regularization,...
            'Shuffle','every-epoch',...
            'Verbose',true,...
            'ExecutionEnvironment',exeEnvironment);
            trainedNet = trainNetwork(valImages,valLabels,trainedNet.Layers,options);
        end

        if ~isempty(valImages)
            predictedLabels = classify(trainedNet,valImages);
            for e = 1:numel(predictedLabels)
                %delta(e) = 1 - corr(double(predictedLabels{e})',double(valLabels{e})');
                delta(e) = norm(double(predictedLabels{e})'-double(valLabels{e})');
            end
            valError = mean(delta);
        else
            valError = 0;

        end
        cons = [];

        
        %{
        fileName = num2str(valError,10) + ".mat";
        save(fileName,'trainedNet','valError','options')
        cons = [];
        %}
end

       
