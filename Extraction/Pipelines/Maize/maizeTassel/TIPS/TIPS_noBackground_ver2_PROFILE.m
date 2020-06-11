function [res] = TIPS_noBackground_ver2( img, oPath, bgFunc ,params )
    res = true;
    
    %TIPS Automated phenotype extraction from tassel images
    %   foreground: string. Path to foreground image
    %   background: string. Path to background image
    %   out: this is now a path only - to match the other functions

    %{
    %   out: string. Path and prefix for output.  File extensions will be added
    %        to this, e.g. out = './test' will yield files ./test_processed.png and 
    %        ./test_out.txt
    %}

    if ~isdeployed()
        addpath(strcat(pwd, '/TIPS_functions'));
        addpath(strcat(pwd, '/background_removal_functions'));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this was blocked off via nate the destroyer
    % all bow down - april, 28 2020
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % created
    mmkdir(oPath);
    [p,nm,ext] = fileparts(img);
    out = [oPath '{originalName_' nm '}'];

    %{
    dirToMake = regexp(out, '(.*)/.*', 'tokens');
    dirToMake = char(dirToMake{:});
    mkdir(dirToMake);
    %}


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this was blocked off via nate the destroyer
    % all bow down - april, 28 2020
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    % If not passed in on the command line, initialize parameter values:
    % params = [gThresh padSize 
    if ~exist('params', 'var')
       params = [ 0.08, 200, 55, 15, 10e-8, 75, 301, .2 ];
    end
    P = nameParams(params);
    reSZ = 10; % resizing parameter - work into param list later...
    boxExpansion = 0.5; % Same with this parameter...
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % created
    if nargin < 4;params = [ 0.08, 200, 55, 15, 10e-8, 75, 301, .2 ];end
    P = nameParams(params);
    reSZ = 10; % resizing parameter - work into param list later...
    boxExpansion = 0.5; % Same with this parameter...
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    %{

           FOR TESTING START HERE:
           img = '/Users/jlg374/Desktop/TIPS_testing/successful_images/DDPSC_02069-2-90_0233.jpg';

           foreground = './601153_rep3_08.54.52_side_green.thumb.jpg';
           background = './601153_rep3_08.54.52_side_green_background.thumb.jpg';
           out = './601153_rep3_08.54.52_side_green';

    %}


    % Initialize variables
    area = NaN;
    BNc = NaN;
    TLadj = NaN;
    tort = NaN;
    compact = NaN;
    fracDim = NaN;
    skelLength = NaN;
    perimeter = NaN;

    try

        fprintf('*****************************************\n')
        fprintf(strcat('Beginning analysis on ', img, '\n'))
        fprintf('*****************************************\n')

        fprintf('\nSYSTEM ARCHITECTURE: ');
        system('uname -m');
        fprintf('\n');


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % this was blocked off via nate the destroyer
        % all bow down - april, 28 2020
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Subtract the background from the image with the tassel in it
        %fprintf('\n\nExtracting tassel and creating Binary\n')
        %tassel = double(imread(img)) / 255;
        % make the binary mask
        %tBin = thresholdTasselImage(tassel, reSZ);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % this was blocked off via nate the destroyer
        % all bow down - april, 28 2020
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [tasselMask,tasselImage] = bgFunc(img,oPath);



        for maskN = 1:size(tasselMask,3)
            tBin = tasselMask(:,:,maskN);


            tBin = fliplr(tBin');
            perimeter = sum(sum(bwmorph(tBin, 'remove')));
            
            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % this was blocked off via nate the destroyer
            % all bow down - april, 28 2020
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Skip if any part of tassel touches the edge of the image
            if max(tBin(1,:)) > 0 || ...
                    max(tBin(:,1)) > 0 || ...
                    max(tBin(size(tBin, 1), :)) > 0 || ...
                    max(tBin(:, size(tBin, 2))) > 0
                msgID = 'cleanBinary:touchesEdge';
                msg = 'binary tassel touches edge of image.  Excluding from analysis';
                exception = MException(msgID, msg);
                throw(exception);
            end
            %}
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % remove edges
            for r = 1:4
                tBin(:,1) = 0;
                tBin = imrotate(tBin,90);
            end


            % Pad image with 0s to ensure indices of downstream analyses don't
            % fall of the edges of the array.
            tBin = padarray(tBin, [P.padSize P.padSize]);
            tassel = padarray(tasselImage, [P.padSize P.padSize]);





            % Smooth and re-threshold
            fprintf('Smoothing...\n')
            tSmooth = smoothTassel(tBin, P.smoothSigma, P.smoothKernelDim);



            % Get convex hull, convex area, and convex image of smoothed tassel
            tSmoothProps = regionprops(tSmooth, 'ConvexHull', 'ConvexArea', 'ConvexImage');

            % Get endpoints and splines
            fprintf('Making endpoints, splines, and all that good stuff...\n')
            [endpoints, branchpoints, splines, ~, spike, skelLength, base] = tasselSkel(tSmooth, P.skelTol, P.skelMinBranch);

            % Identify lowest branch by looking at
            % change in thickness of the base of the spike.
            % Fit spline to and calculate length for Tassel Length
            firstBranch = findSpikeStart( P.spikeWidth, tBin, spike, P.spikeTol);
            if ~isempty(firstBranch)
                firstBranchAdj = find(spike(:,1) == firstBranch, 1, 'last');
                [truncSpline, TLadj] = calcAdjSpikeLength(firstBranchAdj, spike);
            end

            % Calculate tortuosity
            tort = calcTort(truncSpline);

            % Calculate branch number by circle method
            [BNc, circle, lowBranch, radius]  = BNcircle(tBin, spike, branchpoints, [], firstBranchAdj);
            BNc = max(BNc);

            % Calculate fractal dimension using box-counting method
            fprintf('Calculating fractal dimension...\n')
            %fracDim = fractalDim(tBin);
            fracDim = 0;

            % Create image with some analyses visualized
            fprintf('Drawing a picture!\n')
            close all;
            preview = figure('Visible','off');
            plotTassel(tBin, tSmoothProps, splines, endpoints, base, circle, lowBranch, truncSpline);
            set(gca, 'Visible', 'off', 'position', [0 0 1 1], 'units', 'normalized');
            set(preview, 'PaperUnits', 'centimeters');
            set(preview, 'PaperPosition', [0 0 18 12])


            imgOutName = [out '{maskNumber_' num2str(maskN) '}{fileType_processedImage}.png'];
            saveas(preview,imgOutName, 'png');
            close all;



            % Pheno measurements:
            % Get tassel area and compactness
            area = sum(sum(tBin));
            compact = area / tSmoothProps.ConvexArea;

            phenos = {img area BNc TLadj tort compact fracDim skelLength perimeter ''};



            % Write phenotypes
            fprintf('Writing phenos...\n')
            %fileOut = strcat(out, '_out.txt');
            fileOut = [out '{maskNumber_' num2str(maskN) '}{fileType_csvPhenoList}.txt'];
            fprintf(['writing to file:' fileOut '\n']);
            cell2csv(fileOut,phenos);
            %fprintf(['Done with ' foreground '\n'])
            close all;
        end

    catch ME

        % Return error info
        getReport(ME)
        phenos = {img area BNc TLadj tort compact fracDim skelLength perimeter ME.message};
        res = false;
    end
    
end


