function [motion, logger] = estSelfMotion(flow,W,f,opt)
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
%   The idea comes from:
%   Rieger, J. and Lawton, D.T. (1985). Processing Differential Image 
%       Motion, Journal of Optical Society America 2:354–359.
%   We compute the local flow derivative in flow direction to approximate
%   the translational flow. Only those derivative vectors are selected that
%   have a length above threshold calculated as the mean of all lengths.
%   This translational flow field is used to estimate the linear velocity.
%   The rotational velocity is estimated using the bilinear motion
%   constraint given the translational velocity.
%
%   Copyright (C) 2012  Florian Raudies, 04/12/2012, Boston University.
%   License, GNU GPL, free software, without any warranty.
%

% *************************************************************************
% Select a robust estimation method, if specified.
% *************************************************************************
if opt.mfunc
    warning('estSelfMotion:InappropriateInputArgument', ...
        'Option mfunc is not supported!');
end
if opt.ransac
    warning('estSelfMotion:InappropriateInputArgument', ...
        'Option ransac is not supported!'); 
end
if opt.em
    warning('estSelfMotion:InappropriateInputArgument', ...
        'Option em is not supported!');
end
% *************************************************************************
% The esimtation routine starts here.
% *************************************************************************

if size(flow.Dx,1)==1 || size(flow.Dx,2)==1,
    error('estSelfMotion:InappropriateInputArgument',...
          'Flow input needs to be 2D!');
end

% Create new logger.
logger = newEstSelfMotionLogger();
% Fetch data from flow structure.
X   = flow.X;
Y   = flow.Y;
Dx  = flow.Dx;
Dy  = flow.Dy;
% Calculate the derivative in flow direction.
Phi = atan2(Dy, Dx);
Dx_dx = circshift(Dx,[+0 -1]) - Dx;
Dx_dy = circshift(Dx,[-1 +0]) - Dx;
Dy_dx = circshift(Dy,[+0 -1]) - Dy;
Dy_dy = circshift(Dy,[-1 +0]) - Dy;
% Set boundary values to zero.
Dx_dx(:,end) = 0;
Dx_dx(end,:) = 0;
Dx_dy(:,end) = 0;
Dx_dy(end,:) = 0;
Dy_dx(:,end) = 0;
Dy_dx(end,:) = 0;
Dy_dy(:,end) = 0;
Dy_dy(end,:) = 0;
Dx_phi = Dx_dx.*cos(Phi) + Dx_dy.*sin(Phi);
Dy_phi = Dy_dx.*cos(Phi) + Dy_dy.*sin(Phi);
Len = hypot(Dx_phi, Dy_phi);
% Index for derivatives that are larger than the mean computed over the 
% length of all derivatives.
Index = Len>mean(mean(Len));
% Define translational flow as the directional derivatives of the input
% flow.
flowTrans = newImgFlow(X(Index),Y(Index),Dx_phi(Index),Dy_phi(Index));
Omega = [0 0 0];
if length(W)==1, 
    Vel = estTranslationBilinear(flowTrans, W, f, Omega, opt.unbias);
else
    Vel = estTranslationBilinear(flowTrans, W(Index), f, Omega, opt.unbias);
end
Omega = estRotationBilinear(flow, W, f, Vel);
motion = newMotion(Vel, Omega);
