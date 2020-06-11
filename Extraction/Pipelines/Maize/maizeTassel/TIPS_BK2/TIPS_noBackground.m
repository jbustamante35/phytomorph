function [] = TIPS_noBackground( foreground, out, params )
%TIPS Automated phenotype extraction from tassel images
%   foreground: string. Path to foreground image
%   background: string. Path to background image
%   out: string. Path and prefix for output.  File extensions will be added
%        to this, e.g. out = './test' will yield files ./test_processed.png and 
%        ./test_out.txt
if ~isdeployed()
    addpath(strcat(pwd, '/TIPS_functions'));
end

dirToMake = regexp(out, '(.*)/.*', 'tokens');
dirToMake = char(dirToMake{:});
mkdir(dirToMake);

% If not passed in on the command line, initialize parameter values:
% params = [gThresh padSize 
if ~exist('params', 'var')
   params = [ 0.08, 200, 55, 15, 10e-8, 75, 301, .2 ];
end
P = nameParams(params);
reSZ = 10; % resizing parameter - work into param list later...
boxExpansion = 0.5; % Same with this parameter...

%{
       
       FOR TESTING START HERE:
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
    fprintf(strcat('Beginning analysis on ', foreground, '\n'))
    fprintf('*****************************************\n')
    
    fprintf('\nSYSTEM ARCHITECTURE: ');
    system('uname -m');
    fprintf('\n');

    % Subtract the background from the image with the tassel in it
    fprintf('\n\nExtracting tassel and creating Binary\n')
    tassel = double(imread(foreground)) / 255;
    box = findTasselCropBox(tassel, reSZ, boxExpansion);
    boxTassel = imcrop(tassel, box);
    boxMask = thresholdTasselImage(boxTassel, reSZ);
    tBin = tasselMask2OriginalSize(boxMask, box, tassel);
    tBin = fliplr(tBin');
    perimeter = sum(sum(bwmorph(tBin, 'remove')));
    
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
    
    % Pad image with 0s to ensure indices of downstream analyses don't
    % fall of the edges of the array.
    tBin = padarray(tBin, [P.padSize P.padSize]);
    tassel = padarray(tassel, [P.padSize P.padSize]);
    
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
    if ~isempty(firstBranch);
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
    fracDim = fractalDim(tBin);
    
    % Create image with some analyses visualized
    fprintf('Drawing a picture!\n')
    close all;
    preview = figure('Visible','off');
    imgOut = flattenMaskOverlay(flip(permute(tassel,[2 1 3]),2),logical(tBin));
    plotTassel(imgOut, tSmoothProps, splines, endpoints, base, circle, lowBranch, truncSpline);
    set(gca, 'Visible', 'off', 'position', [0 0 1 1], 'units', 'normalized');
    set(preview, 'PaperUnits', 'centimeters');
    set(preview, 'PaperPosition', [0 0 18 12])
    saveas(preview, strcat(out, '_processed.png'), 'png');
    close all;
    
    % Pheno measurements:
    % Get tassel area and compactness
    area = sum(sum(tBin));
    compact = area / tSmoothProps.ConvexArea;
    
    phenos = {foreground area BNc TLadj tort compact fracDim skelLength perimeter ''};
    
catch ME
    % Return error info
    getReport(ME)
    phenos = {foreground area BNc TLadj tort compact fracDim skelLength perimeter ME.message};
end



% Write phenotypes
fprintf('Writing phenos...\n')
fileOut = strcat(out, '_out.txt');



tmpDoc = [];
tmpDoc = generatePhenotypeNode(tmpDoc,foreground,{'imageName'},'imageName');
tmpDoc = generatePhenotypeNode(tmpDoc,area,{'Area'},'Area');
tmpDoc = generatePhenotypeNode(tmpDoc,BNc,{'BranchNumber'},'BranchNumber');
tmpDoc = generatePhenotypeNode(tmpDoc,TLadj,{'TasselLength'},'TasselLength');
tmpDoc = generatePhenotypeNode(tmpDoc,tort,{'Tortuosity'},'Tortuosity');
tmpDoc = generatePhenotypeNode(tmpDoc,compact,{'Compactness'},'Compactness');
tmpDoc = generatePhenotypeNode(tmpDoc,fracDim,{'FractalDimension'},'FractalDimension');
tmpDoc = generatePhenotypeNode(tmpDoc,skelLength,{'SkeletonLength'},'SkeletonLength');
tmpDoc = generatePhenotypeNode(tmpDoc,perimeter,{'PerimeterLength'},'PerimeterLength');
JSON_string = savejson('earDoc',tmpDoc);

%%%%%%%%%%%%%%%%%%%%%%%%%
% save JSON string
%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Writing JSON file. \n');
JSONfile = strrep(fileOut,'.txt','.json');
fileID = fopen(JSONfile,'w');
fprintf(fileID,strrep(JSON_string,'\/','\\/'));
fclose(fileID);


fprintf(['writing to file:' fileOut '\n']);
%{
fileID=fopen(fileOut, 'w');
formatSpec=strcat('''%s''', repmat(['\t%f'], 1, size(phenos, 2)-2),'\t''%s''');
fprintf(fileID, formatSpec, phenos{1,:});
fclose(fileID);
%}

cell2csv(fileOut,phenos);
%fprintf(['Done with ' foreground '\n'])
close all;
end


