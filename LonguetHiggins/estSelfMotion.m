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
%   See the following references:
%   Longuet-Higgins, H.C. (1981). A Computer Algorithm for Reconstructing 
%       a Scene from Two Projections. Nature 293:133-135.
%   Hartley, R. and Zisserman, A. (2003). Multiple View Geometry in 
%       Computer Vision. 2nd Edition. Cambridge University Press 
%       (Chapter 9).
%   Hartley, R. (1997). In Defense of the Eight-Point Algorithm. IEEE 
%       Transactions on Pattern Analysis and Machine Intelligence
%       19(6):580-593.
%   In our this implementation we assume the point correspondeces to be
%   given by (X-Dx/2,Y-Dy/2) and (X+Dx/2,Y+Dy/2). For the estimation of the
%   translational velocity we use the normalization proposed by Hartley.
%   The rotational velocity is estimated using the bilinear constraint.
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
% Compute the two viewpoints using the image flow field.
X1 = W.*(X-Dx/2); 
Y1 = W.*(Y-Dy/2);
X2 = W.*(X+Dx/2);
Y2 = W.*(Y+Dy/2);
% Apply transformation of normalization to data points.
[A1, B1, T1] = normalizeData(X1, Y1, f);
[A2, B2, T2] = normalizeData(X2, Y2, f);
% Compute constraint matrix, where f is set to one explicitly, because of
% the applied normalization.
C = [A1.*A2, B1.*A2, A2, ...
     A1.*B2, B1.*B2, B2, ...
     A1, B1, ones(size(A1))]';
[Vec Val] = eig(C*C');
[~, mi] = min(diag(Val));
% Recover translational vector and rotation matrix.
E = reshape(Vec(:,mi), [3 3])';
% Apply transform (because E is calculated for the inverse transform).
E = T1'*E*T2;
[U, ~, ~] = svd(E);
Vel     = U(:,3)';
% Assume a forward translational velocity.
Vel     = Vel*sign(Vel(3));
Omega   = estRotationBilinear(flow, W, f, Vel);
% Create new motion structure with estimates.
motion  = newMotion(Vel, Omega);


function [X, Y, T] = normalizeData(X, Y, f)
% Compute the transfomration to normalize the data and apply it.
meanX = mean(X(:));
meanY = mean(Y(:));
meanD = sqrt(1/(2*numel(X))*sum( (X(:)-meanX).^2 + (Y(:)-meanY).^2 ));
X = (X - meanX) / meanD;
Y = (Y - meanY) / meanD;
T = [1/meanD  0       -meanX/meanD/f; ...
     0        1/meanD -meanY/meanD/f; ...
     0        0        1/f];
