function [motion, logger, A] = estSelfMotionGivenDepth(flow, W, f, opt)
% estSelfMotionGivenDepth
%   flow    - Structure with fields X, Y, Dx, and Dy (see newImgFlow).
%   W       - Weights that are included into the estimation process.
%   f       - Focal length of the pinhole camera.
%   opt     - Structure with options (see newEstSelfMotionOpt).
%             Note that depth is assumed in the field 'Z'.
%
% RETURN
%   motion  - Structure with fields for the linear velocity "Vel" and 
%             rotational velocity "Omega".
%
% DESCRIPTION
%   See the following references:
%   Perrone, J.A. (1992). Model for the Computation of Self-Motion in 
%       Biological Systems. Optical Society of America A, 9(2):177–192.
%   Perrone, J.A. and Stone, L.S. (1994). A Model of Self-Motion Estimation 
%       within Primate Extrastriate Visual Cortex. Vision Research 
%       34(21):2917–2938.
%   This method applies the voting idea (templates) from Perrone & Stone.
%   In this algorithm we assume the depths to be given which reduces the 
%   voting space from 6D to 5D.
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
type = opt.tuneOptimizeType;
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
if strcmp(type,'hg')
    [motion, logger, A] = voteSelfMotionHG(flow,W,f,opt);    
else
    [motion, logger, A] = voteSelfMotion(flow,W,f,opt);
end


% *************************************************************************
% Private function for the estimation of self-motion for a voting map.
% *************************************************************************
function [motion, logger, A] = voteSelfMotion(flow, W, f, opt)
MapPhi      = opt.MapPhi;
MapRad      = opt.MapRad;
MapOmegaX   = opt.MapOmegaX;
MapOmegaY   = opt.MapOmegaY;
MapOmegaZ   = opt.MapOmegaZ;
sigmaSpeed  = opt.tuneSigmaSpeed;
MapVelX = cos(MapPhi).*MapRad;
MapVelY = sin(MapPhi).*MapRad;
MapVelN = sqrt(MapVelX.^2 + MapVelY.^2 + 1);
MapVelX = MapVelX./MapVelN;
MapVelY = MapVelY./MapVelN;
MapVelZ = 1./MapVelN;
% Get dimensions
vecNum  = numel(flow.X);
phiNum  = size(MapPhi,1);
radNum  = size(MapRad,2); 
oxNum   = size(MapOmegaX,1); 
oyNum   = size(MapOmegaX,2); 
ozNum   = size(MapOmegaX,3);
InvZ    = 1./(opt.Z+eps);
% Replicate image coordinates, flow, self-motions, and depths for 6D space.
X   = repmat(flow.X(:),  [1 phiNum radNum oxNum oyNum ozNum]);
Y   = repmat(flow.Y(:),  [1 phiNum radNum oxNum oyNum ozNum]);
Dx  = repmat(flow.Dx(:), [1 phiNum radNum oxNum oyNum ozNum]);
Dy  = repmat(flow.Dy(:), [1 phiNum radNum oxNum oyNum ozNum]);
InvZ = repmat(InvZ(:),   [1 phiNum radNum oxNum oyNum ozNum]);
VelX = shiftdim(repmat(MapVelX, [1 1 oxNum oyNum ozNum vecNum]), 5);
VelY = shiftdim(repmat(MapVelY, [1 1 oxNum oyNum ozNum vecNum]), 5);
VelZ = shiftdim(repmat(MapVelZ, [1 1 oxNum oyNum ozNum vecNum]), 5);
OmegaX = shiftdim(repmat(MapOmegaX, [1 1 1 vecNum phiNum radNum]), 3);
OmegaY = shiftdim(repmat(MapOmegaY, [1 1 1 vecNum phiNum radNum]), 3);
OmegaZ = shiftdim(repmat(MapOmegaZ, [1 1 1 vecNum phiNum radNum]), 3);
if length(W)>1, W = repmat(W, [1 phiNum radNum oxNum oyNum ozNum]); end
% For computational efficiency organize these 6D spaces as linear data 
% structures.
X       = X(:);
Y       = Y(:);
Dx      = Dx(:);
Dy      = Dy(:);
W       = W(:);
InvZ    = InvZ(:);
VelX    = VelX(:);
VelY    = VelY(:);
VelZ    = VelZ(:);
OmegaX  = OmegaX(:);
OmegaY  = OmegaY(:);
OmegaZ  = OmegaZ(:);
% Define template flow field.
Du = (-VelX*f + VelZ.*X).*InvZ ...
   + (X.*Y.*OmegaX -(f^2+X.^2).*OmegaY +f*Y.*OmegaZ)/f;
Dv = (-VelY*f + VelZ.*Y).*InvZ ...
   + ((f^2+Y.^2).*OmegaX -X.*Y.*OmegaY -f*X.*OmegaZ)/f;
% Calculate the similarity by using a rectified cos-tuning function for the 
% directional difference and a exp-log^2 tuning function for the speed
% difference of input flow and template flow vectors.
FlowLen = hypot(Dx, Dy);
TemLen  = hypot(Du, Dv);
% Compute the inner product between input and template.
A       = (Dx.*Du + Dy.*Dv)./(FlowLen.*TemLen+eps);
A(A<0)  = 0;
A       = A.*exp(-(log2(FlowLen./(TemLen+eps)+eps)).^2 / (2*sigmaSpeed^2));
A       = W.*A;
% Reshape the matrix A to be 7D.
A = reshape(A, [vecNum phiNum radNum oxNum oyNum ozNum]);
% Compute the mean for all flow vectors within the visual field. Note that
% the flow templates are defined globally and assume a rigid world.
A = squeeze(mean(A));
% Calculate the index for the maximum match.
[mv, mi] = max(A(:));
[iPhi iRad iOx iOy iOz] = ind2sub([phiNum radNum oxNum oyNum ozNum], mi);
% Determine the translational and rotational motion based on this index.
Vel     = [MapVelX(iPhi,iRad) MapVelY(iPhi,iRad) MapVelZ(iPhi,iRad)];
Omega   = [MapOmegaX(iOx,iOy,iOz) ...
           MapOmegaY(iOx,iOy,iOz) ...
           MapOmegaZ(iOx,iOy,iOz)];
% Create new motion structure with estimates.
motion = newMotion(Vel, Omega);
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
phiNum  = 7;
radNum  = 7;
oxNum   = 5;
oyNum   = 5;
ozNum   = 5;
% Define the initial heading map with the above parameters.
[MapPhi,MapRad] = ndgrid(linspace(MapPhi(1,1),MapPhi(end,1),phiNum),...
                         linspace(MapRad(1,1),MapRad(1,end),radNum));
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
for iter = 1:iterNum-1,
    L = L0.*[(1/(phiNum*iter)+(1/phiNum)^iter),...
             (1/(radNum*iter)+(1/radNum)^iter),...
             (1/(oxNum*iter)+(1/oxNum)^iter),...
             (1/(oyNum*iter)+(1/oyNum)^iter),...
             (1/(ozNum*iter)+(1/ozNum)^iter)]/2;
    [motion, ~, A] = voteSelfMotion(flow,W,f,opt);
    [~, mi] = max(A(:));
    [iPhi iRad , ~, ~, ~] = ind2sub([phiNum radNum oxNum oyNum ozNum], mi);
    % Define a new voting space that is centered around the maximum.
    [MapPhi,MapRad] = ...
        ndgrid(linspace(MapPhi(iPhi,iRad)-L(1), MapPhi(iPhi,iRad)+L(1), phiNum),...
               linspace(MapRad(iPhi,iRad)-L(2), MapRad(iPhi,iRad)+L(2), radNum));
    [MapOmegaX,MapOmegaY,MapOmegaZ] = ...
        ndgrid(linspace(motion.Omega(1)-L(3), motion.Omega(1)+L(3), oxNum),...
               linspace(motion.Omega(2)-L(4), motion.Omega(2)+L(4), oyNum),...
               linspace(motion.Omega(3)-L(5), motion.Omega(3)+L(5), ozNum));
    % Save the new voting space in the opt(ions) structure.
    opt.MapPhi      = MapPhi;
    opt.MapRad      = MapRad;
    opt.MapOmegaX   = MapOmegaX;
    opt.MapOmegaY   = MapOmegaY;
    opt.MapOmegaZ   = MapOmegaZ;   
end
% One estimation outside the loop to estimate depths as well.
[motion, logger, A] = voteSelfMotion(flow, W, f, opt);
% Book-keeping of iterations.
logger.initMean  = 1;
logger.iterMean  = iterNum;
