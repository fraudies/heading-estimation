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
%   Kanatani, K. (1993). 3-D Interpretation of Optical Flow by 
%       Renormalization. International Journal of Computer Vision, 
%       11(3):276-282.
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

% *************************************************************************
% Setup of auxiliary data structures.
% *************************************************************************
% Calculate m.
M_n = sqrt(X.^2+Y.^2+f^2)+eps;
M1  = flow.X(:)./M_n;
M2  = flow.Y(:)./M_n;
M3  = f./M_n;
% Calculate m star dot.
M_tmp = +(X.*Dx+Y.*Dy)./(M_n.^3);
Mdt1  = +Dx./M_n-M_tmp.*X;
Mdt2  = +Dy./M_n-M_tmp.*Y;
Mdt3  = -f*M_tmp;
% Calculate vector cross product Mst = M x Mdt.
Mst1 = M2.*Mdt3 - M3.*Mdt2;
Mst2 = M3.*Mdt1 - M1.*Mdt3;
Mst3 = M1.*Mdt2 - M2.*Mdt1;
% Calculate L.
l11 = sum(W.*Mst1.^2);
l12 = sum(W.*Mst1.*Mst2);
l13 = sum(W.*Mst1.*Mst3);
l22 = sum(W.*Mst2.^2);
l23 = sum(W.*Mst2.*Mst3);
l33 = sum(W.*Mst3.^2);
L   = [l11 l12 l13;...
       l12 l22 l23;...
       l13 l23 l33];
% Calculate M1ij, M2ij, M3ij.
M11 = W.*M1.^2;
M12 = W.*M1.*M2;
M13 = W.*M1.*M3;
M22 = W.*M2.^2;
M23 = W.*M2.*M3;
M33 = W.*M3.^2;
m111 = sum(Mst1.*M11);
m112 = sum(Mst1.*M12);
m113 = sum(Mst1.*M13);
m122 = sum(Mst1.*M22);
m123 = sum(Mst1.*M23);
m133 = sum(Mst1.*M33);
M16  = [m111, m112, m122, m113, m123, m133];

m211 = sum(Mst2.*M11);
m212 = sum(Mst2.*M12);
m213 = sum(Mst2.*M13);
m222 = sum(Mst2.*M22);
m223 = sum(Mst2.*M23);
m233 = sum(Mst2.*M33);
M26  = [m211, m212, m222, m213, m223, m233];

m311 = sum(Mst3.*M11);
m312 = sum(Mst3.*M12);
m313 = sum(Mst3.*M13);
m322 = sum(Mst3.*M22);
m323 = sum(Mst3.*M23);
m333 = sum(Mst3.*M33);
M36  = [m311, m312, m322, m313, m323, m333];

% Calcualte N, see Appendix D of Kanatani (1993).
N = zeros(6,6);
N(1,1) = sum(W.*M1.^4);
N(1,2) = sum(W.*M1.^3.*M2);
N(1,3) = sum(W.*M1.^2.*M2.^2);
N(1,4) = sum(W.*M1.^3.*M3);
N(1,5) = sum(W.*M1.^2.*M2.*M3);
N(1,6) = sum(W.*M1.^2.*M3.^2);

N(2,1) = N(1,2);
N(2,2) = sum(W.*M1.^2.*M2.^2);
N(2,3) = sum(W.*M1.*M2.^3);
N(2,4) = sum(W.*M1.^2.*M2.*M3);
N(2,5) = sum(W.*M1.*M2.^2.*M3);
N(2,6) = sum(W.*M1.*M2.*M3.^2);

N(3,1) = N(1,3);
N(3,2) = N(2,3);
N(3,3) = sum(W.*M2.^4);
N(3,4) = sum(W.*M2.^2.*M1.*M3);
N(3,5) = sum(W.*M2.^3.*M3);
N(3,6) = sum(W.*M2.^2.*M3.^2);

N(4,1) = N(1,4);
N(4,2) = N(2,4);
N(4,3) = N(3,4);
N(4,4) = sum(W.*M1.^2.*M3.^2);
N(4,5) = sum(W.*M1.*M2.*M3.^2);
N(4,6) = sum(W.*M1.*M3.^3);

N(5,1) = N(1,5);
N(5,2) = N(2,5);
N(5,3) = N(3,5);
N(5,4) = N(4,5);
N(5,5) = sum(W.*M2.^2.*M3.^2);
N(5,6) = sum(W.*M2.*M3.^3);

N(6,1) = N(1,6);
N(6,2) = N(2,6);
N(6,3) = N(3,6);
N(6,4) = N(4,6);
N(6,5) = N(5,6);
N(6,6) = sum(W.*M3.^4);

% Transformation to compute the inverse.
N(:, [2 4 5]) = 2*N(:, [2 4 5]);

% *************************************************************************
% Caluculate the translational velocity.
% *************************************************************************
B16 = N\M16';
B26 = N\M26';
B36 = N\M36';

% From 1 x 6 representation to 1 x 9 representation!
B1 = [B16(1) B16(2) B16(4); B16(2) B16(3) B16(5); B16(4) B16(5) B16(6)];
B2 = [B26(1) B26(2) B26(4); B26(2) B26(3) B26(5); B26(4) B26(5) B26(6)];
B3 = [B36(1) B36(2) B36(4); B36(2) B36(3) B36(5); B36(4) B36(5) B36(6)];

B19 = B1(:);
B29 = B2(:);
B39 = B3(:);

M19 = [m111; m112; m113;  m112; m122; m123;  m113; m123; m133];
M29 = [m211; m212; m213;  m212; m222; m223;  m213; m223; m233];
M39 = [m311; m312; m313;  m312; m322; m323;  m313; m323; m333];

A = L - [M19'*B19 M19'*B29 M19'*B39;...
         M29'*B19 M29'*B29 M29'*B39;...
         M39'*B19 M39'*B29 M39'*B39];
     
% Reduce statistical bias.
if opt.unbias,
    m11 = sum(M11);
    m12 = sum(M12);
    m13 = sum(M13);
    m22 = sum(M22);
    m23 = sum(M23);
    m33 = sum(M33);
    M   = eye(3)*length(flow.X(:)) ...
        - [m11 m12 m13;  m12 m22 m23;  m13 m23 m33];
    Msqrt       = M^(-0.5);
    [Vec, Val]  = eig(Msqrt*A*Msqrt);
    [~, mi]     = min(diag(Val));
    Vel         = Vec(:,mi)';
    Vel         = Vel*Msqrt;
else 
    [Vec, Val]  = eig(A);
    [~, mi]     = min(diag(Val));
    Vel         = Vec(:,mi)';
end
% Assume forward translation.
Vel = Vel/(Vel(3)+eps);
% Assume unit-speed.
Vel = Vel/norm(Vel);

% *************************************************************************
% Calculate rotation.
% *************************************************************************
K       = -(B1*Vel(1)+B2*Vel(2)+B3*Vel(3));
Omega   = ( 0.5 * ( (trace(K)+3*Vel*K*Vel')*Vel' ) - 2*K*Vel' )';

% Create new motion structure with estimates.
motion  = newMotion(Vel, Omega);
