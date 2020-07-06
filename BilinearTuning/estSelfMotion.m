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
%   This method uses the bilinear constraint to estimate the self-motion
%   with a voting method. The maximum log-likelihood estimate for 
%       exp(-(bilinear constraint)^2)
%   is computed.
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

% *************************************************************************
% The esimtation routine starts here.
% *************************************************************************
type            = opt.tuneOptimizeType;
% Map for directions of translational velocity in polar coordinates.
[MapPhi, MapRad] = ndgrid(opt.TunePhi, opt.TuneRad);
[MapOmegaX,MapOmegaY,MapOmegaZ] = ...
    ndgrid(opt.TuneOmegaX, opt.TuneOmegaY, opt.TuneOmegaZ);
% Add these maps to the structure opt.
opt.MapPhi      = MapPhi;
opt.MapRad      = MapRad;
opt.MapOmegaX   = MapOmegaX;
opt.MapOmegaY   = MapOmegaY;
opt.MapOmegaZ   = MapOmegaZ;
% Apply a hierarchical grid method with the voting.
if strcmp(type, 'hg')
    [motion, logger, A] = voteSelfMotionHG(flow,W,f,opt);
else
    [motion, logger, A] = voteSelfMotion(flow,W,f,opt);
end


% *************************************************************************
% Private function for the estimation of self-motion for a given voting 
% map.
% *************************************************************************
function [motion, logger, A] = voteSelfMotion(flow,W,f,opt)
% Fetch the data of the flow.
X   = flow.X(:);
Y   = flow.Y(:);
Dx  = flow.Dx(:);
Dy  = flow.Dy(:);
% Retrieve information about the maps from the opt(ions) structure.
MapPhi      = opt.MapPhi;
MapRad      = opt.MapRad;
MapOmegaX   = opt.MapOmegaX;
MapOmegaY   = opt.MapOmegaY;
MapOmegaZ   = opt.MapOmegaZ;
% Get the dimensions for the voting space.
[phiNum radNum]     = size(MapPhi);
[oxNum oyNum ozNum] = size(MapOmegaX);
vecNum              = numel(X);
% Configure the map of possible translational velocities.
MapVelX = cos(MapPhi).*MapRad;
MapVelY = sin(MapPhi).*MapRad;
% Setup the voting space by creating 5D matrices.
VelX    = shiftdim(repmat(MapVelX,[1 1 oxNum oyNum ozNum vecNum]), 5);
VelY    = shiftdim(repmat(MapVelY,[1 1 oxNum oyNum ozNum vecNum]), 5);
VelZ    = 1;
OmegaX  = shiftdim(repmat(MapOmegaX,[1 1 1 vecNum phiNum radNum]), 3);
OmegaY  = shiftdim(repmat(MapOmegaY,[1 1 1 vecNum phiNum radNum]), 3);
OmegaZ  = shiftdim(repmat(MapOmegaZ,[1 1 1 vecNum phiNum radNum]), 3);
X       = repmat(X, [1 phiNum radNum oxNum oyNum ozNum]);
Y       = repmat(Y, [1 phiNum radNum oxNum oyNum ozNum]);
Dx      = repmat(Dx, [1 phiNum radNum oxNum oyNum ozNum]);
Dy      = repmat(Dy, [1 phiNum radNum oxNum oyNum ozNum]);
% For computational efficiency organize these as linear data structures.
X       = X(:);       
Y       = Y(:);      
Dx      = Dx(:);   
Dy      = Dy(:);
VelX    = VelX(:); 
VelY    = VelY(:); 
OmegaX  = OmegaX(:); 
OmegaY  = OmegaY(:); 
OmegaZ  = OmegaZ(:);
% Use the bilinear constraint with rho = 0 as similarity constraint for the 
% voting.
A = (VelX.*(-f*Dy -X.*Y.*OmegaY +(f^2+Y.^2).*OmegaX -f*X.*OmegaZ) ...
    +VelY.*(+f*Dx -f*Y.*OmegaZ +(f^2+X.^2).*OmegaY -X.*Y.*OmegaX) ...
    +VelZ.*(+X.*Dy -Y.*Dx -f*X.*OmegaX+(X.^2+Y.^2).*OmegaZ -f*Y.*OmegaY)).^2;
% Free up memory.
clear X Y Dx Dy VelX VelY OmegaX OmegaY OmegaZ;
% If we have non-uniform weights apply them.
if length(W)>1
    W = repmat(W, [1 phiNum radNum oxNum oyNum ozNum]);
    W = W(:);
end
% Estimate log-likelihood.
A = -W.*A;
A = reshape(A, [vecNum phiNum radNum oxNum oyNum ozNum]);
A = squeeze(mean(A, 1));
% Read-out activity.
[mv, mi] = max(A(:));
[iPhi iRad iOx iOy iOz] = ind2sub([phiNum radNum oxNum oyNum ozNum], mi);
Vel     = [MapVelX(iPhi,iRad) MapVelX(iPhi,iRad) 1];
Omega   = [MapOmegaX(iOx,iOy,iOz) MapOmegaY(iOx,iOy,iOz) MapOmegaZ(iOx,iOy,iOz)];
% Normalize translational velocity.
Vel     = Vel/norm(Vel);
% Define structure with motion.
motion  = newMotion(Vel,Omega);
% Create new logger.
logger              = newEstSelfMotionLogger();
logger.maxActivity  = mv;
logger.phiNum       = phiNum;
logger.radNum       = radNum;
logger.oxNum        = oxNum;
logger.oyNum        = oyNum;
logger.ozNum        = ozNum;

% *************************************************************************
% Private function for estimation of self-motion for a given map. This
% method only uses the boundaries of the map to apply a hierarchical grid
% (hg) search.
% *************************************************************************
function [motion, logger, A] = voteSelfMotionHG(flow,W,f,opt)
% Retrieve information about the maps from the opt(ions) structure.
MapPhi      = opt.MapPhi;
MapRad      = opt.MapRad;
MapOmegaX   = opt.MapOmegaX;
MapOmegaY   = opt.MapOmegaY;
MapOmegaZ   = opt.MapOmegaZ;
% Set prameters for the adaptive voting spaces.
iterNum = 10;
phiNum  = 5;
radNum  = 5;
oxNum   = 5;
oyNum   = 5;
ozNum   = 5;
% Define the initial heading map with the above parameters.
[MapPhi,MapRad] = ndgrid(linspace(MapPhi(1),MapPhi(end),phiNum),...
                         linspace(MapRad(1),MapRad(end),radNum));
[MapOmegaX,MapOmegaY,MapOmegaZ] = ...
    ndgrid(linspace(MapOmegaX(1,1,1),MapOmegaX(end,1,1),oxNum),...
    	   linspace(MapOmegaY(1,1,1),MapOmegaY(1,end,1),oyNum),...
           linspace(MapOmegaZ(1,1,1),MapOmegaZ(1,1,end),ozNum));
% Define the length of each interval of the 5D space.
L0 = [MapPhi(end,1) - MapPhi(1,1),...
      MapRad(1,end) - MapRad(1,1),...
      MapOmegaX(end,1,1) - MapOmegaX(1,1,1),...
      MapOmegaY(1,end,1) - MapOmegaY(1,1,1),...
      MapOmegaZ(1,1,end) - MapOmegaZ(1,1,1)];
% Store the map parameters in the structure opt.
opt.MapPhi      = MapPhi;
opt.MapRad      = MapRad;
opt.MapOmegaX   = MapOmegaX;
opt.MapOmegaY   = MapOmegaY;
opt.MapOmegaZ   = MapOmegaZ;
% Iterative search by using a hierarchical grid with overlapping intervals
% between iterations.
for iter = 1:iterNum,
    % Define the interval length for this iteration.
    L = L0.*[(1/(phiNum*iter)+(1/phiNum)^iter), ...
             (1/(radNum*iter)+(1/radNum)^iter), ...
             (1/(oxNum*iter)+(1/oxNum)^iter), ...
             (1/(oyNum*iter)+(1/oyNum)^iter), ...
             (1/(ozNum*iter)+(1/ozNum)^iter)]/2;
    % Do the voting using the above method.
    [motion, ~, A] = voteSelfMotion(flow,W,f,opt);
    % Find the propotype with the most votes.
    [~, mi]     = max(A(:));
    [iPhi iRad , ~, ~, ~] = ind2sub([phiNum radNum oxNum oyNum ozNum], mi);
    % Define a new voting space that is centered around the maximum.
    [MapPhi,MapRad] = ...
        ndgrid(linspace(MapPhi(iPhi,iRad)-L(1),MapPhi(iPhi,iRad)+L(1),phiNum),...
               linspace(MapRad(iPhi,iRad)-L(2),MapRad(iPhi,iRad)+L(2),radNum));
    [MapOmegaX,MapOmegaY,MapOmegaZ] = ...
        ndgrid(linspace(motion.Omega(1)-L(3), motion.Omega(1)+L(3), oxNum),...
               linspace(motion.Omega(2)-L(4), motion.Omega(2)+L(4), oyNum),...
               linspace(motion.Omega(3)-L(5), motion.Omega(3)+L(5), ozNum));
    % Save the new voting space in the opt(ions) structure.
    opt.MapPhi = MapPhi;
    opt.MapRad = MapRad;
    opt.MapOmegaX = MapOmegaX;
    opt.MapOmegaY = MapOmegaY;
    opt.MapOmegaZ = MapOmegaZ;
end
% Book-keeping of iterations.
logger.initMean  = 1;
logger.iterMean  = iterNum;
