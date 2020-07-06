function [motion, logger] = estSelfMotion(~, ~, ~, opt)
% estSelfMotion
%   flow    - Structure with fields X, Y, Dx, and Dy (see newImgFlow).
%   W       - Weights that are included into the estimation process.
%   f       - Focal length of the pinhole camera.
%   opt     - Structure with options (see newEstSelfMotionOpt).
%
% RETURN
%   motion  - Structure with fields for the linear velocity "Vel" and 
%             rotational velocity "Omega".
%
% DESCRIPTION
%   A most common self-motion for humans is the locomotion into the
%   direction of gaze and no rotations -- a straight walking. Formally,
%   this is expressed as the self-motion Vel = [0 0 1] and Omega = [0 0 0]
%   which is returned by this method as a prior for heading.
%
%   Copyright (C) 2012  Florian Raudies, 04/12/2012, Boston University.
%   License, GNU GPL, free software, without any warranty.
%

% *************************************************************************
% Select a robust estimation method, if specified.
% *************************************************************************
if opt.ransac, 
    warning('estSelfMotion:InappropriateInputArgument', ...
        'Option ransac is not supported!'); 
end
if opt.mfunc, 
    warning('estSelfMotion:InappropriateInputArgument', ...
        'Option mfunc is not supported!'); 
end
if opt.em,
    warning('estSelfMotion:InappropriateInputArgument', ...
        'Option em is not supported!'); 
end

% Create new logger.
logger       = newEstSelfMotionLogger();
% Create new motion structure with estimates.
motion = newMotion([0 0 1], [0 0 0]);
