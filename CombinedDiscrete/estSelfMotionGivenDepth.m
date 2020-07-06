function [motion, logger] = estSelfMotionGivenDepth(flow, W, f, opt)
% estSelfMotionGivenDepth
%   flow    - Structure with fields X, Y, Dx, and Dy (see newImgFlow).
%   W       - Weights that are included into the estimation process.
%   f       - Focal length of the pinhole camera.
%   opt     - Structure with options (see newEstSelfMotionOpt).
%             Note that depth is assumed in the field 'Z'.
%
% RETURN
%   motion  - Structure with fields for the linear velocity "Vel" and 
%             rotational velocity "Omega".
%
% DESCRIPTION
%   This method estimates self-motion in one step considering depths given.
%
%   Copyright (C) 2012  Florian Raudies, 04/12/2012, Boston University.
%   License, GNU GPL, free software, without any warranty.
%

% *************************************************************************
% Select a robust estimation method, if specified.
% *************************************************************************
if opt.mfunc
    [motion, logger] = estSelfMotionMFunc(flow, f, opt); return;
end
if opt.ransac
    [motion logger] = estSelfMotionRANSAC(flow, f, opt); return;
end
if opt.em
    warning('estSelfMotion:InappropriateInputArgument', ...
        'Option em is not supported!'); 
end

% *************************************************************************
% The esimtation routine starts here.
% *************************************************************************
% Create new logger.
logger = newEstSelfMotionLogger();
if isfield(opt,'Index'), opt.Z = opt.Z(opt.Index); end
InvZ   = 1./(opt.Z+eps);
motion = estSelfMotionPseudoInverse(flow, W, f, InvZ);
