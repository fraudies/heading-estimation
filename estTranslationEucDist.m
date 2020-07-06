function Vel = estTranslationEucDist(flow, W, f, Omega, InvZ)
% estTranslationEucDist
%   flow    - Structure with fields X, Y, Dx, and Dy (see newImgFlow).
%   W       - Weights that are included into the estimation process.
%   f       - Focal length of the pinhole camera.
%   Omega   - Rotational velocity.
%   InvZ    - Inverse of the depth Z.
%
% RETURN
%   Vel     - Linear velocity.
%
% DESCRIPTION
%   Least-squares optimization with squared Euclidean difference between 
%   input flow and model flow:
%        _
%       /_ W |(Dx,Dy) - 1/Z A Vel - B Omega|^2 -Vel-> min
%
%   Copyright (C) 2012  Florian Raudies, 04/12/2012, Boston University.
%   License, GNU GPL, free software, without any warranty.
%

% Load flow field into column-vectors.
X    = flow.X(:);
Y    = flow.Y(:);
Dx   = flow.Dx(:);
Dy   = flow.Dy(:);
InvZ = InvZ(:);
% Estimate the translational velocity knowing the flow, inverse depth, and
% rotational velocity.
A1  = W.*(Dx -1/f*(X.*Y*Omega(1) -(f^2+X.^2)*Omega(2) +f*Y*Omega(3)));
A2  = W.*(Dy -1/f*((f^2+Y.^2)*Omega(1) -X.*Y*Omega(2) -f*X*Omega(3)));
J   = [sum(-f*InvZ.*A1), sum(-f*InvZ.*A2), sum(InvZ.*(X.*A1+Y.*A2))];
e1  = +f^2*sum(W.*InvZ.^2);
e2  = -f*sum(W.*InvZ.^2.*X);
e3  = -f*sum(W.*InvZ.^2.*Y);
e4  = +sum(W.*InvZ.^2.*(X.^2+Y.^2));
E   = [e1  0 e2;...
       0  e1 e3;...
       e2 e3 e4];
Vel = J/E;
% Assume forward translation.
Vel = Vel/(Vel(3)+eps);
% Assume unit-speed.
Vel = Vel/norm(Vel);
