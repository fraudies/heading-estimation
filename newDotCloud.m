function [scene camera] = newDotCloud(camera, zMin, zMax, xNum, yNum)
% newDotCloud
%   camera  - Camera model with inverse projection function and image
%             plane.
%   zMin    - Minmum depth along the optical axis in METER.
%   zMax    - Maximum depth in METER.
%   xNum    - Number of samples in the horizontal dimension.
%   yNum    - Number of samples in the vertical dimension.
%   
% RETURN
%   scene   - Scene defined by 3D coordinates X, Y, and Z.
%   camera  - Camera with updated image plane and zBuffer.
%
% DESCRIPTION
%   This method defines first the image plane with its samples, then the
%   depth coordinates. Second, it uses the inverse projection method of the
%   camera model to convert the image coordinates into world coordinates.
%   This technique allows to define the 3D dot cloud as a cuboid in 
%   projective space. In the 3D world this is 3D dot cloud is a frustum.
%
%   Copyright (C) 2012  Florian Raudies, 04/12/2012, Boston University.
%   License, GNU GPL, free software, without any warranty.
%

% Get the dimensions of the image plane from the camera model.
xMin = -camera.width/2;
xMax = +camera.width/2;
yMin = -camera.height/2;
yMax = +camera.height/2;
% Define the image plane of the camera.
[Y X] = ndgrid(linspace(yMin,yMax,yNum),linspace(xMin,xMax,xNum));
imgPlane = newImgPlane(X,Y);
camera.imgPlane = imgPlane;
% Fill the zBuffer of the camera.
camera.Z = rand(yNum,xNum)*(zMax-zMin)+zMin;
% Use the inverse projection method to define the 3D coordinates.
scene = camera.fInvProj(camera);
