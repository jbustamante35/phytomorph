function [] = renderMap(data,overlayImage,skipVec,displayValues,iVec,mM,oPath)
    try

        % for each waveThreshold
        for waveLevel = 1:size(data,3)
            MIN_LEVEL = mM(waveLevel,1);
            MAX_LEVEL = mM(waveLevel,2);

            % make the heat map of the X
            tmp = data(:,:,waveLevel);

            if iVec(1)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % make the heat map
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % mask overlay
                out = zeros(size(overlayImage));
                % get the color map
                CM = jet(MAX_LEVEL);
                % flip the colormap
                CM = flipdim(CM,1);
                cnt = 1;
                for f = 1:MAX_LEVEL
                    % reshape the inFrame wave transistion
                    overLayMask = round(tmp) == f;
                    out = flattenMaskOverlay(out,overLayMask,.75,CM(cnt,:));
                    %imshow(out,[]);
                    %drawnow
                    cnt = cnt + 1;
                end
                imwrite(out,[oPath '{heatMap_waveThreshold_' num2str(waveLevel) '}.tif']);
            end

        end


         for waveLevel = 1:size(data,3)
            MIN_LEVEL = mM(waveLevel,1);
            MAX_LEVEL = mM(waveLevel,2);

            % make the heat map of the X
            tmp = data(:,:,waveLevel);


            if iVec(2)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % make the t-stack
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % for each t-SKIP
                for skip = 1:numel(skipVec)
                    
                    % set the current skip
                    SKIP = skipVec(skip);
                    % create the level vec
                    levelVec = linspace(MIN_LEVEL,MAX_LEVEL,round((MAX_LEVEL-MIN_LEVEL)/SKIP));
                    % show the first frame
                    imshow(overlayImage,displayValues);hold on
                    % get the number of frames
                    NF = numel(levelVec);
                    % get the color map for NF number of frames
                    CM = jet(NF);

                   



                    % for each level in the level vector
                    for level = 1:numel(levelVec)
                        % get the contours at each of the PER_contour levels
                        dC = contourc(data(:,:,waveLevel),[levelVec(level) levelVec(level)]);
                        % init variables for the contours
                        cnt = 1;CON = {};
                        % while there are still contours
                        while ~isempty(dC)
                            % get the contour pointer
                            ptr = dC(:,1);
                            % get the contour
                            CON{cnt} = dC(:,2:1+ptr(2));
                            % pop/remove the contour from the queue
                            dC(:,1:1+ptr(2)) = [];
                            plot(CON{cnt}(1,:),CON{cnt}(2,:),'Color',CM(level,:));hold on
                            cnt = cnt + 1;
                        end
                        drawnow
                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % save the results of the SKIP t-projections
                    saveas(gca,[oPath '{skipLevel_' num2str(SKIP) '}' '{contourMap_' num2str(waveLevel) '}.tif']);
                    hold off
                end
            end


         end
    catch ME
        ME
    end
    
end