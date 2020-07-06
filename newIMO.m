function imo = newIMO(motion, Pos, Shape)
% newIMO
%   motion  - Motion of the IMO.
%   Pos     - Position (x,y) of the IMO in the image plane in PERCENT.
%   Shape   - Size (w,h) of the IMO in the image plane in PERCENT.
%
% RETURN
%   imo     - Structure with the above parameters as fields.
%
%   Copyright (C) 2012  Florian Raudies, 04/12/2012, Boston University.
%   License, GNU GPL, free software, without any warranty.
%

imo.Pos     = Pos;
imo.Shape   = Shape;
imo.name    = 'IMO';
imo.motion  = motion;
imo.fHandle = @addIMO2Flow;


% *************************************************************************
% Private function of IMO structures.
% *************************************************************************
function flow = addIMO2Flow(camera, noise, flow)
% addIMO2Flow
%   camera  - Information of the camera model.
%   noise   - Information about the indpendently moving object.
%   flow    - Image flow field. 
%
% RETURN
%   flow    - Flow field with motion for IMO replaced in the provided flow.
%
%   Copyright (C) 2012  Florian Raudies, 04/12/2012, Boston University.
%   License, GNU GPL, free software, without any warranty.

if ~strcmp(noise.name,'IMO'), error('Only IMO noise for this method'); end
% Retrieve position and shape of the IMO from the noise structure.
Pos = noise.Pos;
Shape = noise.Shape;
% Fetch image coordiantes for this camera model.
X = camera.imgPlane.X;
Y = camera.imgPlane.Y;
% Fetch zBuffer for this scene.
Z = camera.Z;
xPos    = camera.width*Pos(1);
yPos    = camera.height*Pos(2);
if any(Shape==0), error('Invalid shape in newIMO'); end
width   = camera.width*Shape(1);
height  = camera.height*Shape(2);
IMO     = abs(X+camera.width/2-xPos)<=width/2 ...
        & abs(Y+camera.height/2-yPos)<=height/2;
% Restrict camera to IMO region.
camera.imgPlane.X   = X(IMO);
camera.imgPlane.Y   = Y(IMO);
camera.Z            = Z(IMO);
% Compute the flow for the IMO.
flowIMO         = camera.fImgFlow(camera, noise.motion);
% Replace the original flow by the IMO flow.
flow.X(IMO)     = flowIMO.X;
flow.Y(IMO)     = flowIMO.Y;
flow.Dx(IMO)    = flowIMO.Dx;
flow.Dy(IMO)    = flowIMO.Dy;

