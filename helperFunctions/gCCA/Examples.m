
%% Example:  samples as the common mode

% The data set is the "flavoured waters" data described in:
% * Måge I, Menichelli E, Næs T. Preference mapping by PO-PLS: Separating common and unique information in several data blocks. Food Qual Prefer. 2012;24(1):8-16.
clear all; close all;
load waters

X={Odour Flavour Liking}
options=PCAGCA(X);
options.BlockNames=strvcat('Odour','Flavour','Liking');

% Example 1: Automatic selection of components, based on limits for
% canonical correlation and explained variance
options.ExpVarLim=5; %each component has to explain at least 5% of the variation in the blocks
options.Rlim=0.90; %the canonical correlation has to be above 0.9 in order to define a component as common
options.autoselect=1;
model=PCAGCA(X,options)

% The allocation of scores is given in model.Labels:
model.Labels{1} %shows if the scores are common (C) or unique/distinctive (U)

%scatter plot of scores
figure; sample_grouping(model.Scores{3}(:,1),model.Scores{3}(:,2),Liking.i(:,6),Liking.i)

% Example 2: Predefine the number of principal components for each block,
% but interactivly decide how many common:
options=PCAGCA(X)
options.nCompsLocalPCA=[3 3 6];
model=PCAGCA(X,options)


% Example 3: Interactive selection of all components:
% A number of plots will pop up, that are useful for selecting the numbers of components. 
% In each case, you will have to enter the number of components in the command window to continue. 
% The component selection is sequential:
% 1. Decide a number of principal components for each block (based on separate PCA models for each block)
% 2. For each common subspace, decide the number of components based on canonical correlation coefficient (the dots in the plot) and explained variance in each bloc (bars)
model=PCAGCA(X,options)



%% Example:  variables as the common mode

% The data set is the simulated, with two global common components
clear all; close all
load simdata

X={X1 X2 X3}
options=PCAGCA(X)
options.nCompsLocalPCA=[5 5 5];
model=PCAGCA(X,options)

