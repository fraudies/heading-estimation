function InvZ = estInvZ(flow, f, motion)
% estInvZ
%   flow    - Structure with fields X, Y, Dx, and Dy (see newImgFlow).
%   f       - Focal length of the pinhole camera.
%   motion  - Structure with fields for the linear velocity "Vel" and 
%             rotational velocity "Omega".
%
% RETURN
%   InvZ    - Inverse depth Z estimated from flow and self-motion.
%
% DESCRIPTION
%   Depths are estimated employing a least squares optimization using 
%   the quadratic Euclidean difference between input and model flow. 
%     _
%    /_  |(flw.U(:),flw.V(:)) - 1/Z A Vel - B Omega|^2 -Z-> min
%
%   Copyright (C) 2012  Florian Raudies, 04/12/2012, Boston University.
%   License, GNU GPL, free software, without any warranty.
%

% Fetch motion parameters.
Vel     = motion.Vel;
Omega   = motion.Omega;
% Fetch flow parameters in row-order.
X   = flow.X(:);
Y   = flow.Y(:);
Dx  = flow.Dx(:);
Dy  = flow.Dy(:);
% Compute the residual after subtraction of the rotational flow.
DxRes = Dx - 1/f*(X.*Y*Omega(1) -(f^2+X.^2)*Omega(2) +f*Y*Omega(3));
DyRes = Dy - 1/f*((f^2+Y.^2)*Omega(1) -X.*Y*Omega(2) -f*X*Omega(3));
% Compute the inverse depth.
InvZ  = (DxRes.^2 + DyRes.^2) ...
     ./ (-f*DxRes*Vel(1) -f*DyRes*Vel(2) +(X.*DxRes+Y.*DyRes)*Vel(3) + eps);
