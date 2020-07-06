function Vel = estTranslationBilinear(flow, W, f, Omega, unbias)
% estTranslationBilinear
%   flow    - Structure with fields X, Y, Dx, and Dy (see newImgFlow).
%   W       - Weights that are included into the estimation process.
%   f       - Focal length of the pinhole camera.
%   Omega   - Rotational velocity.
%   unbias  - Flag that indicates if the statistical bias should be
%             reduced.
%
% RETURN
%   Vel     - Linear velocity.
%
% DESCRIPTION
%   Least-squares optimization with the bilinear constraint:
%     _
%    /_ W [(A Vel)^t ( (Dx,Dy) - B Omega)]^2 -Vel-> min
%
%   This bilinear constraint estimation is biased by |AT|^2.
%
%   Copyright (C) 2012  Florian Raudies, 04/12/2012, Boston University.
%   License, GNU GPL, free software, without any warranty.
%

% Fetch data from flow structure.
X   = flow.X(:);
Y   = flow.Y(:);
Dx  = flow.Dx(:);
Dy  = flow.Dy(:);
% Caluclate translation.
C1 = -f*Dy  -X.*Y*Omega(2) +(f^2+Y.^2)*Omega(1) -f*X*Omega(3);
C2 = +f*Dx  -f*Y*Omega(3)  +(f^2+X.^2)*Omega(2) -X.*Y*Omega(1);
C3 = +X.*Dy -Y.*Dx -f*X*Omega(1) +(X.^2+Y.^2)*Omega(3) -f*Y*Omega(2);
c11 = sum(W.*C1.^2);
c22 = sum(W.*C2.^2);
c33 = sum(W.*C3.^2);
c12 = sum(W.*C1.*C2);
c13 = sum(W.*C1.*C3);
c23 = sum(W.*C2.*C3);
C = [c11 c12 c13;...
     c12 c22 c23;...
     c13 c23 c33];
% Check for bias correction.
if unbias,
    g1 = +f^2*mean(W.^2);
    g2 = -f*mean(W.^2.*X);
    g3 = -f*mean(W.^2.*Y);
    g4 = +mean(W.^2.*(X.^2+Y.^2));
    G = [g1 0  g2; ...
         0  g1 g3; ...
         g2 g3 g4];
    G = G^(-1/2);
    [Vec, Val]  = eig(G*C*G);
    [~, i]      = min(diag(Val));
    Vel         = G*Vec(:,i);
else    
    [Vec, Val]  = eig(C);
    [~, i]      = min(diag(Val));
    Vel         = Vec(:,i);
end
% Assume translational motion into forward direction.
Vel = Vel/(Vel(3)+eps);
% Assume unit-speed.
Vel = Vel'/norm(Vel);
