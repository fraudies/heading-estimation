function [motion logger] = estSelfMotionRANSAC(flow, f, opt)
% estSelfMotionRANSAC
%   flow    - Image flow field.
%   f       - Focal length of the pinhole camera.
%   opt     - Parameters for the algorithm.
%
% RETURN
%   motion  - Estimated motion (the most probable one).
%   logger  - Data log of overall probability for models and all estimated
%             motions.
%
% DESCRIPTION
%   This method uses the Random Sample Consensus algorithm to estimate a
%   robust model.
%   Fischler, M. and Bolles, R. (1981). Random Sample Consensus: A Paradigm 
%       for Model Fitting with Applications to Image Analysis and Automated 
%       Cartography, Commun. ACM 24(6):381–395.
%
%   Copyright (C) 2012  Florian Raudies, 04/12/2012, Boston University.
%   License, GNU GPL, free software, without any warranty.

% Do not apply any other robust optimization.
opt.em      = 0;
opt.mest    = 0;
opt.ransac  = 0;
% Create new logger.
logger  = newEstSelfMotionLogger();
% Load flow field in columnar organization.
X       = flow.X(:);   
Y       = flow.Y(:);   
Dx      = flow.Dx(:);   
Dy      = flow.Dy(:);
vecNum  = length(X);
% Load parameters of the RANSAC algorithm.
initNum         = opt.ransacInitNum;
triesNum        = opt.ransacTriesNum;
etaVel          = opt.ransacEtaVel;
etaOmega        = opt.ransacEtaOmega;
setNumFrac      = opt.ransacSetNumFrac;
setNum          = ceil(setNumFrac * vecNum + initNum);
% Setup a data structure to hold all estimated models.
Models = zeros(triesNum, 7); % #data, Vel, Omega
NumTestModels = zeros(triesNum, 1);
% All tries.
for k = 1:triesNum
    % Generate a randomized index.
    P = randperm(vecNum);
    % Construct initial index using the first initNum entries of P.
    S = P(1:initNum);
    opt.Index = S; % Add index for selection of corresponding depths.
    motion = opt.fHandle(newImgFlow(X(S),Y(S),Dx(S),Dy(S)), 1, f, opt);
    % Subsequently add new data points to the initial index if these
    % data points are "compatible" with the model.
    for i = initNum+1:vecNum
        S_ = [S P(i)];
        opt.Index = S_; % Add index for selection of corresponding depths.
        motion_ = opt.fHandle(newImgFlow(X(S_),Y(S_),Dx(S_),Dy(S_)), 1 , f, opt);
        NumTestModels(k) = i-initNum;
        if (norm(motion_.Vel-motion.Vel)<etaVel) ...
        && (norm(motion_.Omega-motion.Omega)<etaOmega)
            % Update set and motion model.
            S = S_; motion = motion_;
            if length(S) > setNum, 
                logger.modelType = 'good';
                logger.dataNum  = length(S);
                logger.initMean = k;
                logger.iterMean = mean(NumTestModels(1:k));
                logger.iterStd  = std(NumTestModels(1:k));
                return; 
            end
        end
    end
    Models(k,:) = [length(S) motion.Vel motion.Omega];
end
% Even though no good model was found, return the best model that is the
% model compatible with most data points.
[mv, mi] = max(Models(:,1));
motion = newMotion(Models(mi,2:4), Models(mi,5:7)); % Vel, Omega.
% Book-keeping about iterations and model quaility.
logger.model = 'best of tested';
logger.dataNum  = mv;
logger.initMean = triesNum;
logger.iterMean = mean(NumTestModels);
logger.iterStd  = std(NumTestModels);
