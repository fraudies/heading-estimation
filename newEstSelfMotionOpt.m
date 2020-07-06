function opt = newEstSelfMotionOpt(fHandle, name, label)
% newEstSelfMotionOpt
%   fhandle             - Function handle to algorithm. SEE function.
%   name                - Name of algorithm. Only for the purpose to save.
%   label               - Label to be displayed in legend to identify this
%                         algorithm.
%
% RETURNS
%   opt                - Structure with the above listed fields.
%
% DESCRIPTION
%   Creates structure for identification of algorithm and input/output 
%   arguments when calling this algorithm.
%
%   Copyright (C) 2012  Florian Raudies, 04/12/2012, Boston University.
%   License, GNU GPL, free software, without any warranty.
%

% General information.
opt.fHandle            = fHandle;
opt.name               = name;
opt.label              = label;
opt.fResidual          = @residualBilinear;
opt.unbias             = 1;
opt.mfunc              = 0;
opt.ransac             = 0;
opt.em                 = 0;
% Hierarchical grid optimization (HG) and definition of intervals for 
% voting spaces.
opt.tuneOptimizeType   = 'hg';
opt.TunePhi            = linspace(2*pi/7, 2*pi, 7);
opt.TuneRad            = exp2space(-4, 2, 7);
opt.TuneOmegaX         = linspace(-10, +10, 5)/180*pi;
opt.TuneOmegaY         = linspace(-10, +10, 5)/180*pi;
opt.TuneOmegaZ         = linspace(-10, +10, 5)/180*pi;
opt.TuneDepth          = exp2space(1, 5, 5);
opt.tuneSigmaSpeed     = sqrt(0.5);
% Parameters for fix-point iteration (FP).
opt.fixPtMaxIter       = 500;
opt.fixPtInitNum       = 15;
opt.fixPtEtaVel        = 5*10^-3;
opt.fixPtEtaOmega      = 10^-4;
% Parameters for Gauss Newton iteration (GN).
opt.gnMaxIter          = 500;
opt.gnInitNum          = 15;
opt.gnEtaDVel          = 10^-7;
opt.gnLambdaRho        = 0.25;
% Parameters for optimization with m-function(s).
opt.mFuncMaxIter       = 50;
opt.mFuncEtaW          = 10^-4;
opt.mFuncSigma         = 10^-6;
opt.mFuncRho           = 0;
% Parameters for Random Sample Consensus (RANSAC).
opt.ransacInitNum      = 9;
opt.ransacTriesNum     = 10^2;
opt.ransacEtaVel       = 10^-3;
opt.ransacEtaOmega     = 10^-3;
opt.ransacSetNumFrac   = 1/3;
% Parameters for Expectation Maximization algorithm (EM).
opt.emMaxIter          = 25;
opt.emModelNum         = 2;
opt.emModelSigma       = 10^-2;
opt.emEtaDrop          = 10^-7;
opt.emEtaVel           = 10^-3;
opt.emEtaOmega         = 10^-3;
opt.emEtaPi            = 10^-4;
opt.emEtaSigma         = 10^-4;
