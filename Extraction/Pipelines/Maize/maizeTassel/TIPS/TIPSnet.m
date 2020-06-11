function [Y,Xf,Af] = TIPSnet(X,~,~)
%TIPSNET neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 04-May-2020 11:14:11.
% 
% [Y] = TIPSnet(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 3xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 2xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [14;12;0];
x1_step1.gain = [0.00836820083682008;0.00862068965517241;0.0078740157480315];
x1_step1.ymin = -1;

% Layer 1
b1 = [4.1322837835149739405;-0.36414299696159369013;-3.6036181843540977532;-1.5032219452245016988;-2.3963782552832357808];
IW1_1 = [-5.7988084073158923104 -3.3811495854161632835 2.9678905701146964802;-1.4843711395117704654 -0.0060935271729361225115 2.423998581412078579;-7.0933598991776749543 -8.4215384656227794125 16.01763478592837231;-0.80933280395536999485 -1.6953422266954414344 -0.19945061966840621492;-1.2634583888971502308 -1.2894860786043669254 -1.2758574810682385969];

% Layer 2
b2 = [0.5684973485506354951;-0.9264107412564439592];
LW2_1 = [-1.783324616916708294 6.6139476581923135612 -3.5426217048185368874 0.23638466338603159045 0.04296748758798919765;0.59919482901201492897 -5.9974016098545313014 3.2364483774382732939 -1.0499129161719975567 -0.12190863562801611264];

% ===== SIMULATION ========

% Format Input Arguments
isCellX = iscell(X);
if ~isCellX
  X = {X};
end

% Dimensions
TS = size(X,2); % timesteps
if ~isempty(X)
  Q = size(X{1},2); % samples/series
else
  Q = 0;
end

% Allocate Outputs
Y = cell(1,TS);

% Time loop
for ts=1:TS

    % Input 1
    Xp1 = mapminmax_apply(X{1,ts},x1_step1);
    
    % Layer 1
    a1 = tansig_apply(repmat(b1,1,Q) + IW1_1*Xp1);
    
    % Layer 2
    a2 = softmax_apply(repmat(b2,1,Q) + LW2_1*a1);
    
    % Output 1
    Y{1,ts} = a2;
end

% Final Delay States
Xf = cell(1,0);
Af = cell(2,0);

% Format Output Arguments
if ~isCellX
  Y = cell2mat(Y);
end
end

% ===== MODULE FUNCTIONS ========

% Map Minimum and Maximum Input Processing Function
function y = mapminmax_apply(x,settings)
  y = bsxfun(@minus,x,settings.xoffset);
  y = bsxfun(@times,y,settings.gain);
  y = bsxfun(@plus,y,settings.ymin);
end

% Competitive Soft Transfer Function
function a = softmax_apply(n,~)
  if isa(n,'gpuArray')
    a = iSoftmaxApplyGPU(n);
  else
    a = iSoftmaxApplyCPU(n);
  end
end
function a = iSoftmaxApplyCPU(n)
  nmax = max(n,[],1);
  n = bsxfun(@minus,n,nmax);
  numerator = exp(n);
  denominator = sum(numerator,1); 
  denominator(denominator == 0) = 1;
  a = bsxfun(@rdivide,numerator,denominator);
end
function a = iSoftmaxApplyGPU(n)
  nmax = max(n,[],1);
  numerator = arrayfun(@iSoftmaxApplyGPUHelper1,n,nmax);
  denominator = sum(numerator,1);
  a = arrayfun(@iSoftmaxApplyGPUHelper2,numerator,denominator);
end
function numerator = iSoftmaxApplyGPUHelper1(n,nmax)
  numerator = exp(n - nmax);
end
function a = iSoftmaxApplyGPUHelper2(numerator,denominator)
  if (denominator == 0)
    a = numerator;
  else
    a = numerator ./ denominator;
  end
end

% Sigmoid Symmetric Transfer Function
function a = tansig_apply(n,~)
  a = 2 ./ (1 + exp(-2*n)) - 1;
end