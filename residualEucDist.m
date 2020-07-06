function R2 = residualEucDist(flow, f, motion, InvZ)
% residualEucDist
%   flow    - Image flow field.
%   f       - Focal length of the pinhole camera.
%   motion  - Motion of the camera.
%   InvZ    - Inverse depth.
%
% RETURN
%   R       - Residual value between flow and model for motion.
%
% DESCRIPTION
%   Computes the residual between flow and model using the Euclidian 
%   squared distance:
%
%       R = |(Dx,Dy) - InvZ A Vel - B Omega|^2 
%
%   Copyright (C) 2012  Florian Raudies, 04/12/2012, Boston University.
%   License, GNU GPL, free software, without any warranty.
%

% Fetch motion.
Vel     = motion.Vel;
Omega   = motion.Omega;
% Fetch flow.
X       = flow.X(:);
Y       = flow.Y(:);
Dx      = flow.Dx(:);
Dy      = flow.Dy(:);
InvZ    = InvZ(:);
% Compute the model flow.
DxModel = (-Vel(1)*f + Vel(3).*X).*InvZ ...
        + (X.*Y.*Omega(1) - (f^2+X.^2).*Omega(2) + f*Y.*Omega(3))/f;
DyModel = (-Vel(2)*f + Vel(3).*Y).*InvZ ...
        + ((f^2+Y.^2).*Omega(1) - X.*Y.*Omega(2) - f*X.*Omega(3))/f;
R2 = (DxModel-Dx).^2 +(DyModel-Dy).^2;
