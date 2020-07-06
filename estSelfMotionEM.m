function [motion, logger] = estSelfMotionEM(flow, f, opt)
% estSelfMotionEM
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
%   This method uses the Expectation Maximization (EM) algorithm with 
%   Gaussian Mixture models to estimate the self-motions that are contained
%   in the flow. Per default the method assumes two models. References for
%   further reading are:
%   MacLean, W. Jepson, A. and Frecker, R. (1994). Recovery of Egomotion 
%       and Segmentation of Independent Object Motion using the EM 
%       Algorithm, British Machine Vision Conference 175–184.
%   Dempster, A.P., Laird, N.M., and Rubin, D.B. (1977). Maximum Likelihood 
%       from Incomplete Data via the EM Algorithm. Journal of the Royal 
%       Statistical Society. Series B (Methodological) 39(1):1-38.
%
%   Copyright (C) 2012  Florian Raudies, 04/12/2012, Boston University.
%   License, GNU GPL, free software, without any warranty.
%

% Do not combine with any other robust estimation method.
opt.em          = 0;
opt.mfunc       = 0;
opt.ransac      = 0;
% Create new logger.
logger          = newEstSelfMotionLogger();
% Fetch paramter values for the EM algorithm.
vecNum          = numel(flow.X);
maxIter         = opt.emMaxIter;
modelNum        = opt.emModelNum;
modelSigma      = opt.emModelSigma;
etaDrop         = opt.emEtaDrop;
etaVel          = opt.emEtaVel;
etaOmega        = opt.emEtaOmega;
etaPi           = opt.emEtaPi;
etaSigma        = opt.emEtaSigma;
% Initialize data structures.
%                 pi, sigma, vel, omega, invDepth
paramNum        = 1 + 1 + 3 + 3 + vecNum;
ModelParam      = zeros(modelNum, paramNum);
ModelParam(:,1) = 1;
ModelParam(:,2) = modelSigma;
% Initialize one model.
motion = opt.fHandle(flow, 1, f, opt);
ModelParam(1,3:5) = motion.Vel;
ModelParam(1,6:8) = motion.Omega;
ModelParam(1,8+(1:vecNum)) = estInvZ(flow,f,motion);
ModelPiIter     = inf([maxIter, modelNum]);
% Iterate E-step and M-step.
for iter = 1:maxIter,
    % Store model parameters.
    ModelParamOld = ModelParam;
    
    % *********************************************************************
    % DROP MODEL if its assigned probability is below etaDrop.
    % *********************************************************************
    Index       = find(ModelParam(:,1)>etaDrop);
    ModelParam  = ModelParam(Index,:);
    modelCount  = length(Index);
    ModelRSq    = zeros(modelCount, vecNum);
    ModelPdf    = zeros(modelCount, vecNum);
    ModelPi     = repmat(ModelParam(:,1),[1 vecNum]);
    ModelSigma  = ModelParam(:,2);
    
    % *********************************************************************
    % EXPECTATION STEP
    % *********************************************************************
    for im = 1:modelNum,
        motion  = newMotion(ModelParam(im,3:5),ModelParam(im,6:8));
        InvZ    = ModelParam(im,8+(1:vecNum));
        R2      = residualEucDist(flow, f, motion, InvZ);
        ModelRSq(im,:) = R2;
        ModelPdf(im,:) = 1/(sqrt(2*pi)*ModelSigma(im)) ...
                       * exp( -R2/ModelSigma(im)^2 );
    end
    ProbData = ModelPdf.*ModelPi ...
        ./ (repmat(sum(ModelPdf.*ModelPi,1),[modelCount 1]) + eps);
    ModelParam(:,1) = mean(ProbData, 2); % Update model pi.
    ModelPiIter(iter,1:modelCount) = ModelParam(:,1);
    
    % *********************************************************************
    % MAXIMIZATION STEP
    % *********************************************************************
    ModelParam(:,2) = sqrt( sum(ProbData.*ModelRSq, 2)./(sum(ProbData, 2) + eps) );
    % Update parameters with new probabilities.
    for im = 1:modelCount,
        motion = opt.fHandle(flow, ProbData(im,:)', f, opt);
        ModelParam(im,3:5) = motion.Vel;
        ModelParam(im,6:8) = motion.Omega;
        ModelParam(im,8+(1:vecNum)) = estInvZ(flow,f,motion);
    end
    
    % *********************************************************************
    % INITIALIZATION OF NEW MODELS
    % *********************************************************************    
    if modelCount < modelNum,
        modelCount = modelCount+1;
        ProbData(modelCount,:) = 1-sum(ProbData,1);
        motion = opt.fHandle(flow, ProbData(im,:)', f, opt);
        ModelParamTmp = ModelParam;
        ModelParam = zeros(modelCount,paramNum);
        ModelParam(1:modelCount-1,:) = ModelParamTmp;
        ModelParam(modelCount,1) = 1;
        ModelParam(modelCount,2) = modelSigma;
        ModelParam(modelCount,3:5) = motion.Vel;
        ModelParam(modelCount,6:8) = motion.Omega;
        ModelParam(modelCount,8+(1:vecNum)) = estInvZ(flow,f,motion);
    elseif modelCount == modelNum ...
       && sum(abs(ModelParam(:,1)-ModelParamOld(:,1))) < etaPi ...
       && sum(abs(ModelParam(:,2)-ModelParamOld(:,2))) < etaSigma ...
       && sum(sqrt(sum((ModelParam(:,3:5)-ModelParamOld(:,3:5)).^2,2))) ...
            < etaVel*modelCount ...
       && sum(sqrt(sum((ModelParam(:,6:8)-ModelParamOld(:,6:8)).^2,2))) ...
            < etaOmega*modelCount,
        disp(iter); break;
    end
    
end
% Return the model with highest probability.
[~, mi] = max(ModelParam(:,1));
motion = newMotion(ModelParam(mi,3:5),ModelParam(mi,6:8));
logger.ModelPiIter = ModelPiIter;
% Add all estimated motion models to the logger.
logger.Motion = ModelParam(:,3:8);
