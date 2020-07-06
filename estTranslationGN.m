function [Vel, dVel] = estTranslationGN(flow, W, f, motion, InvZ, rho)
% estTranslationGN
%   flow    - Structure with fields X, Y, Dx, and Dy (see newImgFlow).
%   W       - Weights that are included into the estimation process.
%   f       - Focal length of the pinhole camera.
%   motion  - Linear and rotational velocity. Note that both velocities are
%             estimates e.g. from prior iterations.
%   InvZ    - Inverse of the depth Z.
%   rho     - Regularization parameter that is adapted during iterations.
%
% RETURN
%   Vel     - Modified translation with Gauss-Newton update.
%   dVel    - Length for translational vector of the Gauss-Newton update.
%
% DESCRIPTION
%   Use Gauss-Newton iteration to calculate a translational update. The 
%   schema is defined as [d/dx f(x)]x = x0 detlaX = -f(x0). Here:
%
%   [d/dVel ((Dx,Dy)^t - 1/Z A Vel - B Omega') A'Vel' 
%       / |A'Vel'|^rho ]Vel=Vel' deltaT 
%           = -((Dx,Dy)^t - B Omega')^t(A'Vel')/|A'Vel'|^rho
%
%   with ' denoting the transpose operator.
%
%   Copyright (C) 2012  Florian Raudies, 04/12/2012, Boston University.
%   License, GNU GPL, free software, without any warranty.
%

EPS_INV = 10^4*eps;
Vel     = motion.Vel;
Omega   = motion.Omega;
InvZ    = InvZ(:);
X   = flow.X(:);
Y   = flow.Y(:);
Dx  = flow.Dx(:);
Dy  = flow.Dy(:);
W   = W(:);
% Estimate translational velocity with adaptive "unbiasing" rho.
E1 = (+f^2*Vel(2) -f*Y*Vel(3)) .* InvZ;
E2 = (-f^2*Vel(1) +f*X*Vel(3)) .* InvZ;
E3 = (-f*X*Vel(2) +f*Y*Vel(1)) .* InvZ;
e11 = sum(W.*E1.^2);
e12 = sum(W.*E1.*E2);
e13 = sum(W.*E1.*E3);
e22 = sum(W.*E2.^2);
e23 = sum(W.*E2.*E3);
e33 = sum(W.*E3.^2);
E = [e11 e12 e13;...
     e12 e22 e23;...
     e13 e23 e33];
% Check for stability of inverse.
if rcond(E+EPS_INV) < EPS_INV, dVel=0; return; end
K = W.*((Dx-1/f*(X.*Y*Omega(1)-(f^2+X.^2)*Omega(2)+f*Y*Omega(3))).*(-f*Vel(2)+Y*Vel(3))...
       +(Dy-1/f*((f^2+Y.^2)*Omega(1)-X.*Y*Omega(2)-f*X*Omega(3))).*(f*Vel(1)-X*Vel(3)));
I = [sum(K.*E1) sum(K.*E2) sum(K.*E3)];
H = (-f*Vel(2)+Y*Vel(3)).^2+(f*Vel(1)-X*Vel(3)).^2;
J = sum(H.^(rho/2))*Vel;
IE      = I/(E+EPS_INV);
JE      = J/(E+EPS_INV);
DVel    = IE-(IE*Vel')/(JE*Vel')*JE;
Vel     = Vel+DVel;
dVel    = norm(DVel);
