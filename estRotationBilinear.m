function Omega = estRotationBilinear(flow, W, f, Vel)
% estRotationBilinear
%   flow    - Structure with fields X, Y, Dx, and Dy (see newImgFlow).
%   W       - Weights that are included into the estimation process.
%   f       - Focal length of the pinhole camera.
%   Vel     - Linear velocity.
%
% RETURN
%   Omega   - Rotational velocity.
%
% DESCRIPTION
%   Least-squares optimization for the bilinear constraint:
%        _
%       /_ W [(A Vel)^t ( (Dx,Dy) - B Omega)]^2 -Omega-> min
%
%   yields the solution for omega.
%
%   Copyright (C) 2012  Florian Raudies, 04/12/2012, Boston University.
%   License, GNU GPL, free software, without any warranty.
%

% Fetch data from flow structure.
X   = flow.X(:);
Y   = flow.Y(:);
Dx  = flow.Dx(:);
Dy  = flow.Dy(:);
W   = W(:);
% Calculate rotational velocity from flow and given translation.
E1  = f*X*Vel(3)  -(f^2+Y.^2)*Vel(1)  +X.*Y*Vel(2);
E2  = X.*Y*Vel(1) -(f^2+X.^2)*Vel(2)  +f*Y*Vel(3);
E3  = f*Y*Vel(2)  -(Y.^2+X.^2)*Vel(3) +f*X*Vel(1);
e11 = sum(W.*E1.^2);
e22 = sum(W.*E2.^2);
e33 = sum(W.*E3.^2);
e12 = sum(W.*E1.*E2);
e13 = sum(W.*E1.*E3);
e23 = sum(W.*E2.*E3);
E   = [e11 e12 e13;...
       e12 e22 e23;...
       e13 e23 e33];
C  = W.*(f*(Dx*Vel(2) -Dy*Vel(1)) +Vel(3)*(X.*Dy-Y.*Dx));
d1 = sum(C.*E1);
d2 = sum(C.*E2);
d3 = sum(C.*E3);
Omega = [d1 d2 d3]/E;
