function [motion, logger] = estSelfMotionMFunc(flow, f, opt)
% estSelfMotionMFunc
%   flow    - Image flow field.
%   f       - Focal length of the pinhole camera.
%   opt     - Options for the algorithm.
%
% RETURN
%   motion  - Estimated camera motion.
%   logger  - Logged number of iterations.
%
% DESCRIPTION
%   Uses a (convex) penalty function for the optimization.
%
%   Copyright (C) 2012  Florian Raudies, 04/12/2012, Boston University.
%   License, GNU GPL, free software, without any warranty.

% Do not combine with any other method for robust estimation.
opt.em       = 0;
opt.mfunc    = 0;
opt.ransac   = 0;
% Create new logger.
logger       = newEstSelfMotionLogger();
% Initialize variables.
mFuncMaxIter = opt.mFuncMaxIter;
mFuncEtaW    = opt.mFuncEtaW;
W            = ones(size(flow.X(:)));
iter         = mFuncMaxIter;
invSigmaW2   = 1 / (2 * opt.mFuncSigma^2);
rho          = opt.mFuncRho;
% Iteration with adaptation of weigths.
for it = 1:mFuncMaxIter
    [motion, logger] = opt.fHandle(flow, W, f, opt);
    Wold    = W; 
    Res     = opt.fResidual(flow, f, motion);
    W       = exp(-Res * invSigmaW2) + rho;
    if mean((W(:)-Wold(:)).^2) < mFuncEtaW, iter = it; break; end
end
% Store the number of iterations into the logger.
logger.mFuncIter = iter;
