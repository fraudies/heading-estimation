function camera = newPinholeCamera(height, width, vFov)
% newPinholeCamera
%   height  - Height of the image plane.
%   width   - Width of the image plane.
%   vFov    - Vertical field of view in RAD.
%
% RETURN
%   camera  - Structure that describes a pinhole camera. Note that this
%             structure also incldues the functionality to project, inverse 
%             project and to compute the image flow for a pinhole camera.
%
%   Copyright (C) 2012  Florian Raudies, 04/12/2012, Boston University.
%   License, GNU GPL, free software, without any warranty.
%

camera.height   = height;
camera.width    = width;
camera.vFov     = vFov;
camera.hFov     = vFov * width/height;
camera.f        = height/(2*tan(vFov/2));
camera.fProj    = @perspProj;
camera.fInvProj = @invPerspProj;
camera.fImgFlow = @pinholeImgFlow;
camera.Z        = 0; % zBuffer, used by other methods.
camera.imgPlane = 0; % structure with image plane, used by other methods.




% *************************************************************************
% Private function for perspective projection of a pinhole camera.
% *************************************************************************
function imgPlane = perspProj(scene, camera)
% perspProj
%   scene       - Scene with X, Y, and Z values.
%   camera      - Pinhole camera with focal length f.
%
% RETURN
%   imgPlane    - Image coordinates in the image plane.
imgPlane.X = camera.f*scene.X./scene.Z;
imgPlane.Y = camera.f*scene.Y./scene.Z;





% *************************************************************************
% Private function for inverse perspective projection of a pinhole camera.
% *************************************************************************
function scene = invPerspProj(camera)
% perspProj
%   camera      - Pinhole camera with X, Y, coordiantes of sensor, 
%                 focal length f, and zBuffer Z.
%
% RETURN
%   scene       - Scene with 3D coordinates in X, Y, and Z.
scene.X = (camera.imgPlane.X).*(camera.Z)/camera.f;
scene.Y = (camera.imgPlane.Y).*(camera.Z)/camera.f;
scene.Z = camera.Z;



% *************************************************************************
% Private function to compute the image flow for a pinhole camera.
% *************************************************************************
function flow = pinholeImgFlow(camera, motion)
% pinholeImgFlow
%   camera  - Structure of the pinhole camera which includes the image 
%             plane and zBuffer.
%   motion  - Self-motion of the camera.
%
% RETURNS
%   flow    - Returns an image flow field for this camera and motion
%             according to the instantaneous motion model proposed by 
%             Longuet-Higgins & Prazdny.
%
% H.C. Longuet-Higgins and K. Prazdny (1980). The Interpretation of a 
% Moving Retinal Image. Proceedings of the Royal Society of London, 
% Series B, Biology Sciences 208, 385-397.
%

% Read image coordinates.
X       = camera.imgPlane.X;
Y       = camera.imgPlane.Y;
% Read corresponding depth values from the zBuffer of the camera.
Z       = camera.Z;
f       = camera.f;
% Read the instantaneous motion parameters for the camera.
Vel     = motion.Vel;
Omega   = motion.Omega;
vx      = Vel(1);
vy      = Vel(2);
vz      = Vel(3);
ox      = Omega(1);
oy      = Omega(2);
oz      = Omega(3);
% Normalize the translational velocity.
vn = sqrt(vx^2 + vy^2 + vz^2) + eps;
vx = vx/vn;     vy = vy/vn;     vz = vz/vn;
% Set fields for output image flow.
flow.X  = X;
flow.Y  = Y;
% Instantaneous motion model from Longuet-Higgins & Prazdny (1980).
flow.Dx = (-vx*f + vz*X)./(Z+eps) ...
        + (X.*Y*ox -(f^2+X.^2)*oy +f*Y*oz)/f;
flow.Dy = (-vy*f + vz*Y)./(Z+eps) ...
        + ((f^2+Y.^2)*ox -X.*Y*oy -f*X*oz)/f;
    
    