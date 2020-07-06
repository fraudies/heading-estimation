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
%   See the following reference:
%   Pauwels, K. and Van Hulle, M. (2007). Optimal Instantaneous Rigid 
%       Motion Estimation Insensitive to Local Minima. Computer Vision and
%       Image Understanding, 104(1):77-86.
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
% Load parameters.
gnMaxIter       = opt.gnMaxIter;
gnInitNum       = opt.gnInitNum;
gnEtaDVel       = opt.gnEtaDVel;
lambdaRho       = opt.gnLambdaRho;
LOG10_EPS_RHO   = -13;
% Initialize map with translational velocity to avoid local minima.
[MapVelX MapVelY MapVelZ] = ...
    sph2cart(rand(gnInitNum,1)*2*pi, rand(gnInitNum,1)*pi, 1);
MapVelX(1)  = 0;
MapVelY(1)  = 0;
MapVelZ(1)  = 1;
resMin      = 10^7;
VelMin      = [0 0 0];
OmegaMin    = [0 0 0];
IterBreak   = zeros(gnInitNum, 1);
% Outer iteration that uses different initial values for translational
% velocities.
for ii = 1:gnInitNum,
    rho = 0;
    % Initialize translational velocity with value from the map.
    Vel = [MapVelX(ii) MapVelY(ii) MapVelZ(ii)];
    for it= 1:gnMaxIter,
        % Estimate the rotational velocity.
        Omega   = estRotationBilinear(flow, W, f, Vel);
        % Estimate inverse depth.
        InvZ    = estInvZ(flow, f, newMotion(Vel,Omega));        
        % Estimate translation and update with GN update.
        [Vel, dVel] = estTranslationGN(flow, W, f, newMotion(Vel,Omega), InvZ, rho);
        % Check for convergence of the Gauss-Newton iteration.
        if dVel < gnEtaDVel, IterBreak(ii) = it; break; end
        % Adaption of regularization parameter.
        rho = min(1, rho+lambdaRho*max(0, log10(dVel)/LOG10_EPS_RHO));
    end
    % Assume forward translation.
    Vel = Vel/(Vel(3)+eps);
    % Assume unit-speed for translational velocity.
    Vel = Vel/norm(Vel);
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
% Book-keeping about iterative method.
logger.iterMean  = mean(IterBreak);
logger.iterStd   = std(IterBreak);
logger.initMean  = gnInitNum;
% Create new motion structure with estimates.
motion = newMotion(Vel, Omega);
