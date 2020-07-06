function [motion, logger, A] = estSelfMotion(flow, W, f, opt)
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
%   Lappe, M. and Rauschecker, J.P. (1993). A Neural Network for the 
%       Processing of Optic Flow from Ego-Motion in Man and Higher Mammals. 
%       Neural Computation, 5:374-391.
%   Lappe, M. (1998). A Model of the Combination of Optic Flow and 
%       Extraretinal Eye Movement Signals in Primate Extrastriate Visual 
%       Cortex - Neural Model of Self-Motion from Optic Flow and 
%       Extraretinal Cues. Neural Networks 11:397-414.
%   Note that unlike these two references and the original proposal of the
%   algorithm we do not use a population encoding of the input flow in this
%   implementation. Instead the input flow is given in a vector
%   representation in Cartesian coordinates.
%
%   Copyright (C) 2012  Florian Raudies, 04/12/2012, Boston University.
%   License, GNU GPL, free software, without any warranty.
%

% *************************************************************************
% Select a robust estimation method, if specified.
% *************************************************************************
if opt.ransac, 
    warning('estSelfMotion:InappropriateInputArgument', ...
        'Option ransac is not supported!'); 
end
if opt.mfunc, 
    warning('estSelfMotion:InappropriateInputArgument', ...
        'Option mfunc is not supported!'); 
end
if opt.em,
    warning('estSelfMotion:InappropriateInputArgument', ...
        'Option em is not supported!'); 
end


type = opt.tuneOptimizeType;
% Map for directions of translational velocity in polar coordinates.
[MapPhi, MapRad] = ndgrid(opt.TunePhi, opt.TuneRad);
opt.MapPhi = MapPhi;
opt.MapRad = MapRad;
% Apply a hierarchical grid method with the voting.
if strcmp(type, 'hg')
    warning('estSelfMotion:InappropriateInputArgument', ...
        'Option hg is not supported!'); 
    %[Vel, logger, A] = voteVelHG(flow,W,f,opt);
    [Vel, logger, A] = voteVel(flow,W,f,opt);
else
    [Vel, logger, A] = voteVel(flow,W,f,opt);
end
% Estimate the rotational velocity.
Omega   = estRotationBilinear(flow, W, f, Vel);
% Create new motion structure with estimates.
motion  = newMotion(Vel, Omega);

% *************************************************************************
% Private function for the estimation of the translation velocity for a 
% heading map.
% *************************************************************************
function [Vel, logger, A] = voteVel(flow, W, f, opt)
% Retrieve information about the maps from the opt(ions) structure.
MapPhi = opt.MapPhi;
MapRad = opt.MapRad;
% Fetch the data of the flow.
X   = flow.X(:);
Y   = flow.Y(:);
Dx  = flow.Dx(:);
Dy  = flow.Dy(:);
vecNum = length(X);
% Set up the translational directions to vote for.
VelX = cos(MapPhi).*MapRad;
VelY = sin(MapPhi).*MapRad;
VelZ = ones(size(VelX));
VelN = sqrt(VelX.^2 + VelY.^2 + VelZ.^2) + eps;
% Normalize the translational velocity to unit-speed.
VelX = VelX./VelN;
VelY = VelY./VelN;
VelZ = VelZ./VelN;
phiNum = size(VelX, 1);
radNum = size(VelY, 2);
% Calculate the similarity for all translational velocities using the
% sub-space constraint from
% Heeger, D.J., Jepson, A.D., (1992). Subspace Methods for Recovering
%   Rigid Motion I: Algorithm and Implementation. International Journal of
%   Computer Vision, 95-117.
A = zeros(phiNum, radNum);
C = zeros(2*vecNum, vecNum+3);
F = zeros(1, 2*vecNum);
F(1:2:end) = Dx;
F(2:2:end) = Dy;
for iPhi = 1:phiNum,
    for iRad = 1:radNum,
        vx = VelX(iPhi,iRad);
        vy = VelY(iPhi,iRad);
        vz = VelZ(iPhi,iRad);
        C(1:2:end,:) = [diag(-vx*f + vz*X), X.*Y/f,    -(f+Y.^2/f), Y];
        C(2:2:end,:) = [diag(-vy*f + vz*Y), (f+Y.^2/f),-X.*Y/f    ,-X];
        A(iPhi,iRad) = sum( (W.*(F*ortc(C))).^2 );
    end
end
[mv, mi] = min(A(:));
[iPhi iRad] = ind2sub([phiNum radNum], mi);
Vel = [VelX(iPhi,iRad) VelY(iPhi,iRad) VelZ(iPhi,iRad)];
% Create new logger.
logger              = newEstSelfMotionLogger();
logger.maxActivity  = mv;
logger.phiNum       = phiNum;
logger.radNum       = radNum;

% % *************************************************************************
% % Private function for estimation of self-motion for a heading map using a
% % hierarchical grid method.
% % Note that this approach fails, due to the non-smoothness of the
% % similarity function. Therefore, it has been disabled.
% % *************************************************************************
% function [Vel, logger, A] = voteVelHG(flow, W, f, opt)
% % Set prameters for the adaptive voting spaces.
% iterNum = 10;
% phiNum  = 5;
% radNum  = 5;
% % Define the initial heading map with the above parameters.
% MapPhi = opt.MapPhi;
% MapRad = opt.MapRad;
% [MapPhi,MapRad] = ndgrid(linspace(MapPhi(1),MapPhi(end),phiNum),...
%                          linspace(MapRad(1),MapRad(end),radNum));
% opt.MapPhi = MapPhi;
% opt.MapRad = MapRad;
% % Define the length of each interval of the 2D space.
% L0 = [MapPhi(end,1)-MapPhi(1,1), MapRad(1,end)-MapRad(1,1)];
% for iter = 1:iterNum,
%     % Define the interval length for this iteration.
%     L = L0.*[(1/(phiNum*iter)+(1/phiNum)^iter),...
%              (1/(radNum*iter)+(1/radNum)^iter)]/2;
%     % Do the voting using the above method.
%     [Vel, ~, A] = voteVel(flow, W, f, opt);
%     % Find the propotype with the most votes.
%     [~, mi] = max(A(:));
%     [iPhi iRad] = ind2sub([phiNum radNum], mi);
%     % Define a new voting space that is centered around the maximum.
%     [MapPhi, MapRad] = ndgrid(...
%         linspace(MapPhi(iPhi,iRad)-L(1), MapPhi(iPhi,iRad)+L(1), phiNum),...
%         linspace(MapRad(iPhi,iRad)-L(2), MapRad(iPhi,iRad)+L(2), radNum));
%     opt.MapPhi = MapPhi;
%     opt.MapRad = MapRad;    
% end
% % Create new logger.
% logger           = newEstSelfMotionLogger();
% % Book-keeping of iterations.
% logger.initMean  = 1;
% logger.iterMean  = iterNum;
