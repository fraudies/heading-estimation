function [motion, logger] = estSelfMotion(flow, W, f, opt)
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
%   This method iterates between the estimation of inverse depths and the
%   estimation of self-motion (translational and rotational velocity). 
%   Self-motion is estimated in one step considering depths given.
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
    [motion, logger] = estSelfMotionEM(flow, f, opt); return;
end

% *************************************************************************
% The esimtation routine starts here.
% *************************************************************************
% Create new logger.
logger       = newEstSelfMotionLogger();
% Load prameters about fix-point iteration from opt(ions).
fixPtMaxIter    = opt.fixPtMaxIter;
fixPtInitNum    = opt.fixPtInitNum;
fixPtEtaVel     = opt.fixPtEtaVel;
fixPtEtaOmega   = opt.fixPtEtaOmega;
% Initialize map with translational velocity to avoid local minima.
[MapVelX MapVelY MapVelZ] = ...
    sph2cart(rand(fixPtInitNum,1)*2*pi, rand(fixPtInitNum,1)*pi, 1);
MapVelX(1)  = 0;
MapVelY(1)  = 0;
MapVelZ(1)  = 1;
resMin      = 10^7;
VelMin      = [0 0 0];
OmegaMin    = [0 0 0];
IterBreak   = zeros(fixPtInitNum, 1);
% Outer iteration that uses different initial values for translational
% velocities.
for ii = 1:fixPtInitNum,
    VelOld      = [0 0 0];    
    OmegaOld    = [0 0 0];
    % Initialize translational velocity with value from the map.
    Vel     = [MapVelX(ii) MapVelY(ii) MapVelZ(ii)];
    Omega   = [0 0 0];
    for it= 1:fixPtMaxIter,
        % Estimate inverse depth.
        InvZ = estInvZ(flow, f, newMotion(Vel,Omega));        
        % Estmate translational and rotational velocity.
        motion = estSelfMotionPseudoInverse(flow, W, f, InvZ);
        % Assume forward translation.
        Vel = motion.Vel/(motion.Vel(3)+eps);
        % Assume unit-speed.
        Vel = Vel/norm(motion.Vel);
        % Retrieve the rotational velocity.
        Omega = motion.Omega;
        % Use L1-norm to check for convergence.
        if max(abs(VelOld-Vel)) < fixPtEtaVel ...
        && max(abs(OmegaOld-Omega)) < fixPtEtaOmega,
            IterBreak(ii) = it;
            break;
        end
        % Update linear and rotational velocity.
        VelOld      = Vel; 
        OmegaOld    = Omega;
    end
    % Calculate the residual between the estimated model and the flow data.
    res = mean( residualBilinear(flow, f, newMotion(Vel, Omega)) );
    % Check if this is a better model than previously estimated ones.
    if (res < resMin),
        VelMin      = Vel;
        OmegaMin    = Omega;
        resMin      = res;
    end
end
Vel     = VelMin;
Omega   = OmegaMin;
% Data-log.
logger.iterMean  = mean(IterBreak);
logger.iterStd   = std(IterBreak);
logger.initMean  = fixPtInitNum;
% Create new motion structure with estimates.
motion = newMotion(Vel, Omega);
