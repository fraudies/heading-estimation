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
%   Heeger, D.J. and Jepson, A.D. (1992). Linear Subspace Methods for
%       Recovering Translational Direction. University of Toronto, 
%       Department of Computer Science, Technical Report RBCV-TR-92-40.
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
% Pull flow data from the structure.
X   = flow.X(:);
Y   = flow.Y(:);
Dx  = flow.Dx(:);
Dy  = flow.Dy(:);
W   = W(:);
% Define the polynomial in X and Y and create a subspace.
vecNum = length(X);
F = [ones(vecNum,1) X Y X.*Y X.^2 Y.^2];
[~, ~, Vs] = svd(F');
% Free up some memory.
clear F;
% The rotation-independent constraints are contained in the (K-6) subspace.
C   = Vs(:,7:vecNum)'; % (vecNum-6) x vecNum
Tau = [f*C*(W.*Dy) -f*C*(W.*Dx) C*(W.*(Y.*Dx-X.*Dy))]; % (vecNum-6) x 3
% Subtract the mean value.
Tau = Tau - repmat(mean(Tau),[vecNum-6, 1]);
% Reduce statistical bias.
if opt.unbias
    g1 = +f^2 * mean(W.^2);
    g2 = -f   * mean(W.^2.*X);
    g3 = -f   * mean(W.^2.*Y);
    g4 = mean(W.^2.*(X.^2+Y.^2));
    G = [g1 0  g2;...
         0  g1 g3;...
         g2 g3 g4];
    G = G^(-1/2);
    [Vec Val]   = eig(G*(Tau'*Tau)*G);
    [~, mi]     = min(diag(Val));
    Vel         = G*Vec(:,mi);
else
    % The product of Tau'*Tau is 3 x 3 matrix.
    [Vec Val]   = eig(Tau'*Tau);
    [~, mi]     = min(diag(Val));
    Vel         = Vec(:,mi);
end
% Assume forward translation.
Vel = Vel/(Vel(3)+eps);
% Assume unit-speed.
Vel = Vel'/norm(Vel);
% Calculate rotation.
Omega   = estRotationBilinear(flow, W, f, Vel);
% Create new motion structure with estimates.
motion  = newMotion(Vel, Omega);
