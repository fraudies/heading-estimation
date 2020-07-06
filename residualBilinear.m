function R = residualBilinear(flow, f, motion)
% residualBilinear
%   flow    - Image flow field.
%   f       - Focal length of the pinhole camera.
%   motion  - Motion of the camera.
%
% RETURN
%   R       - Residual value between flow and model for motion.
%
% DESCRIPTION
%   Computes the residual between flow and model using the bilinear
%   constraint:
%
%       R = |(A Vel)^t ((Dx,Dy) - B Omega)|^2
%
%   Copyright (C) 2012  Florian Raudies, 04/12/2012, Boston University.
%   License, GNU GPL, free software, without any warranty.
%

% Fetch motion.
Vel     = motion.Vel;
Omega   = motion.Omega;
% Fetch flow data.
X   = flow.X(:);
Y   = flow.Y(:);
Dx  = flow.Dx(:);
Dy  = flow.Dy(:);
% Compute residual.
C1  = -f*Dy -X.*Y*Omega(2) +(f^2+Y.^2)*Omega(1) -f*X*Omega(3);
C2  = +f*Dx -f*Y*Omega(3)  +(f^2+X.^2)*Omega(2) -X.*Y*Omega(1);
C3  = +X.*Dy -Y.*Dx -f*X*Omega(1) +(X.^2+Y.^2)*Omega(3) -f*Y*Omega(2);
R   = (Vel(1)*C1 +Vel(2)*C2 +Vel(3)*C3).^2;
