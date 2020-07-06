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
%   Raudies, F. and Neumann, H. (2009). An Efficient Linear Method for the 
%       Estimation of Ego-motion from Optical Flow. In J. Denzler, G. 
%       Notni, and H. Süße (Eds.): DAGM 2009, LNCS 5748, pp. 11-20.
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
% Fetch data from flow structure.
W   = W(:);
X   = flow.X(:);
Y   = flow.Y(:);
Dx  = flow.Dx(:);
Dy  = flow.Dy(:);
% *************************************************************************
% (i) Estimate translational motion.
% *************************************************************************
E   = [-W.*(f^2+Y.^2), W.*X.*Y, f*W.*X, -W.*(f^2+X.^2), f*W.*Y, ...
       -W.*(X.^2+Y.^2)];
M1  = +f*W.*Dy;
M2  = -f*W.*Dx;
M3  = +W.*(Y.*Dx - X.*Dy);
C   = E'*E;
H1  = E'*M1;
H2  = E'*M2;
H3  = E'*M3;
L   = [M1-E*(C\H1) M2-E*(C\H2) M3-E*(C\H3)];
% Reduce the statistical bias.
if opt.unbias,
    g1 = +f^2*mean(W.^2);
    g2 = -f*mean(W.^2.*X);
    g3 = -f*mean(W.^2.*Y);
    g4 = +mean(W.^2.*(X.^2+Y.^2));
    G = [g1 0  g2;...
         0  g1 g3;...
         g2 g3 g4];
    G = G^(-1/2);
    [Vec Val] = eig(G * (L'*L) * G);
    [~, mi] = min(diag(Val));
    Vel = G * Vec(:,mi);
else
    [Vec Val] = eig(L'*L);
    [~, mi] = min(diag(Val));
    Vel = Vec(:,mi);
end
% Assume forward motion and unit speed.
Vel = Vel' / Vel(3);
Vel = Vel / norm(Vel);

% *************************************************************************
% (ii) Estimate rotational motion.
% *************************************************************************
J = Vel(1)*H1 + Vel(2)*H2 + Vel(3)*H3;
I = -[J(1) J(2) J(3);...
      J(2) J(4) J(5);...
      J(3) J(5) J(6)]*Vel';
TmulT = [Vel(1)^2 Vel(2)^2 Vel(3)^2 Vel(1)*Vel(2) Vel(1)*Vel(3) Vel(2)*Vel(3)]';
f1 = [C(1)  C(8)  C(15) 2*C(2)      2*C(3)      2*C(9)     ]*TmulT;
f2 = [C(2)  C(10) C(17) C(4)+C(8)   C(5)+C(9)   C(11)+C(16)]*TmulT;
f3 = [C(3)  C(11) C(18) C(5)+C(9)   C(6)+C(15)  C(12)+C(17)]*TmulT;
f4 = [C(8)  C(22) C(29) 2*C(10)     2*C(11)     2*C(23)    ]*TmulT;
f5 = [C(9)  C(23) C(30) C(11)+C(16) C(12)+C(17) C(24)+C(29)]*TmulT;
f6 = [C(15) C(29) C(36) 2*C(17)     2*C(18)     2*C(30)    ]*TmulT;
F = [f1 f2 f3;...
     f2 f4 f5;...
     f3 f5 f6];
Omega = (F\I)';
% Create new motion structure with estimates.
motion = newMotion(Vel, Omega);
