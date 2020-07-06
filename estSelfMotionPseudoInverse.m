function motion = estSelfMotionPseudoInverse(flow, W, f, InvZ)
% estSelfMotionPseudoInverse
%   flow    - Structure with fields X, Y, Dx, and Dy (see newImgFlow).
%   W       - Weights that are included into the estimation process.
%   f       - Focal length of the pinhole camera.
%   InvZ    - Inverse of the depth Z.
%   opt     - Structure with options (see newEstSelfMotionOpt).
%
% RETURN
%   motion  - Structure with fields for the linear velocity "Vel" and 
%             rotational velocity "Omega".
%
% DESCRIPTION
%   The estimation of the self-motion for given depth is an
%   over-constrained problem. 
%
%     (Dx,Dy) = InvZ*(A * Vel + B * Omega) or
%
%     Dx = (InvZ*A(1,1), InvZ*A(1,2), InvZ*A(1,3), B(1,1), B(1,2), B(1,3))
%          * (Vel, Omega)
%     Dy = (InvZ*A(2,1), InvZ*A(2,2), InvZ*A(2,3), B(2,1); B(2,2); B(2,3))
%          * (Vel, Omega)
%
%   We solve this over-constrained problem using a least-squares approach.
%
%   Copyright (C) 2012  Florian Raudies, 04/12/2012, Boston University.
%   License, GNU GPL, free software, without any warranty.
%


% Fetch data from flow structure.
X       = flow.X(:);
Y       = flow.Y(:);
Dx      = flow.Dx(:);
Dy      = flow.Dy(:);
W       = W(:);
InvZ    = InvZ(:);
N       = zeros(length(X), 1);
% All constraints form an over-determined linear equation system A x = B.
A = [-f*InvZ.*W, N, X.*InvZ.*W,       X.*Y.*W/f, -(f^2+X.^2).*W/f, Y.*W; ...
     N, -f*InvZ.*W, Y.*InvZ.*W, (f^2+Y.^2).*W/f,       -X.*Y.*W/f, -X.*W];
B = [Dx.*W; ...
     Dy.*W];
% We solve the system by taking a least-squares solution (pseudo-inverse).
X = pinv(A)*B;
% Do not apply forward motion and unit-speed constraint to Vel, because 
% for given depths there is no ambiguity.
Vel     = X(1:3)'; 
Omega   = X(4:6)';
motion  = newMotion(Vel,Omega);
