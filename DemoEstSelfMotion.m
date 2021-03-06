clc
close all
clear all

% This method privdes a testbed for the estimation of self-motion from
% image flow. The scene is a random dot could. The self-motion is
% parameterized. 
% An independently moving object (IMO) can be added by the flag "addIMO". 
% Gaussian noise can be added by the flag "addNoise".
% The flow can be displayed setting the flag "showFlow" to one.
% Methods for robustness estimation can be activated by setting the flag
% "ransac", "mfunc", or "em".
% The parameter algoIndex selects one of the 13 algorithms. Chose an index 
% between 1 and 13.
%
%   Copyright (C) 2012  Florian Raudies, 04/12/2012, Boston University.
%   License, GNU GPL, free software, without any warranty.
%

algoIndex   = 1;
addIMO      = 1;
addNoise    = 0;
showFlow    = 1;
ransac      = 1;
mfunc       = 0;
% Note that the option EM is unstable for some algorithms and 
% configurations, e.g. if the IMO does not provide enough data points to 
% estimate a model.
em          = 0;

% *************************************************************************
% Setup the camera and scene and calculate image flow.
% *************************************************************************
% Parameters of the pinhole camera.
width   = 10^-2;        % Width of the image plane in meter.
height  = 10^-2;        % Height of the image plane in meter.
vFov    = 30/180*pi;    % Vertical field of view in RAD.
camera  = newPinholeCamera(width, height, vFov);
% Paramters of the random dot cloud.
zMin    = 2;    % Minimum distance in meter.
zMax    = 10;   % Maximum distance in meter.
xNum    = 10;   % Samples in horiztonal direction.
yNum    = 10;   % Samples in vertical direction.
[scene camera] = newDotCloud(camera, zMin,zMax, xNum,yNum);
% Parameters of instantaneous camera motion.
Vel     = [0 0 1];          % Linear velocity in meter/frame.
Omega   = [5 0 0]/180*pi;   % Rotational velocity in RAD/frame.
motion  = newMotion(Vel, Omega);
% Compute the image flow for this camera model.
flow    = camera.fImgFlow(camera, motion);

% *************************************************************************
% Add independently moving object (IMO) to the flow field.
% *************************************************************************
if addIMO,
    VelIMO      = [0 1 0];
    OmegaIMO    = [0 0 0]/180*pi;
    motionIMO   = newMotion(VelIMO, OmegaIMO);
    Pos         = [0.4 0.6];
    Shape       = [0.5 0.5];
    noise       = newIMO(motionIMO,Pos,Shape);
    flow        = noise.fHandle(camera, noise, flow);
end

% *************************************************************************
% Add Gaussian noise to the flow field.
% *************************************************************************
if addNoise,
    noise   = newGaussianNoise(0.01*camera.f, 0.0);
    flow    = noise.fHandle(noise, flow);
end

% *************************************************************************
% Define all algorithms.
% *************************************************************************
algoNum     = 13;
Algo        = cell(algoNum,1);
% 1st algorithm.
algoName    = 'BilinearPolynom';
addpath(['./',algoName,'/']); 
Algo{1}     = newEstSelfMotionOpt(...
                    @estSelfMotion, algoName, 'bilinear polynom');
rmpath(['./',algoName,'/']);
% 2nd algorithm.
algoName    = 'BilinearTuning'; 
addpath(['./',algoName,'/']); 
Algo{2}     = newEstSelfMotionOpt(...
                    @estSelfMotion, algoName, 'bilinear tuning');
rmpath(['./',algoName,'/']);
% 3rd algorithm.
algoName    = 'BrussHorn'; 
addpath(['./',algoName,'/']); 
Algo{3}     = newEstSelfMotionOpt(...
                    @estSelfMotion, algoName, 'Bruss & Horn, 1983');
rmpath(['./',algoName,'/']);
% 4th algorithm.
algoName    = 'ChiusoEtAl';
addpath(['./',algoName,'/']);
Algo{4}     = newEstSelfMotionOpt(...
                    @estSelfMotion, algoName, 'Chiuso et al., 2000');
rmpath(['./',algoName,'/']);
% 5th algorithm.
algoName    = 'CombinedDiscrete';
addpath(['./',algoName,'/']);
Algo{5}     = newEstSelfMotionOpt(...
                @estSelfMotion, algoName, 'combined discrete');
rmpath(['./',algoName,'/']);
% 6th algorithm.
algoName    = 'HeegerJepson';
addpath(['./',algoName,'/']);
Algo{6}     = newEstSelfMotionOpt(...
                @estSelfMotion, algoName, 'Heeger & Jepson, 1992');
rmpath(['./',algoName,'/']);
% 7th algorithm.
algoName    = 'Kanatani';
addpath(['./',algoName,'/']);
Algo{7}     = newEstSelfMotionOpt(...
                @estSelfMotion, algoName, 'Kanatani, 1993');
rmpath(['./',algoName,'/']);
% 8th algorithm.
algoName    = 'LappeRauschecker';
addpath(['./',algoName,'/']);
Algo{8}     = newEstSelfMotionOpt(...
                @estSelfMotion, algoName, 'Lappe & Rauschecker, 1993');
Algo{8}.TunePhi = linspace(2*pi/11, 2*pi, 11);
Algo{8}.TuneRad = exp2space(-7,     4,    11);
rmpath(['./',algoName,'/']);
% 9th algorithm.
algoName    = 'LonguetHiggins';
addpath(['./',algoName,'/']);
Algo{9}     = newEstSelfMotionOpt(...
                @estSelfMotion, algoName, 'Longuet-Higgins, 1981');
rmpath(['./',algoName,'/']);
% 10th algorithm.
algoName    = 'PauwelsVanHulle';
addpath(['./',algoName,'/']);
Algo{10}     = newEstSelfMotionOpt(...
                @estSelfMotion, algoName, 'Pauwels & Van Hulle, 2007');
rmpath(['./',algoName,'/']);
% 11th algorithm.
algoName     = 'Perrone';
addpath(['./',algoName,'/']);
Algo{11}     = newEstSelfMotionOpt(...
                @estSelfMotion, algoName, 'Perrone, 1992');
rmpath(['./',algoName,'/']);
% 12th algorithm.
algoName     = 'RiegerLawton';
addpath(['./',algoName,'/']);
Algo{12}     = newEstSelfMotionOpt(...
                @estSelfMotion, algoName, 'Rieger & Lawton, 1985');
rmpath(['./',algoName,'/']);
% 13th algorithm.
algoName     = 'ZhangTomasi';
addpath(['./',algoName,'/']);
Algo{13}     = newEstSelfMotionOpt(...
                @estSelfMotion, algoName, 'Zhang & Tomasi, 2002');
rmpath(['./',algoName,'/']);


% *************************************************************************
% Test the selected algorithm.
% *************************************************************************
Algo{algoIndex}.ransac  = ransac;
Algo{algoIndex}.mfunc   = mfunc;
Algo{algoIndex}.em      = em; 
fprintf('Running algorithm: %s.\n\n', Algo{algoIndex}.label);
[motionEst logger]      = Algo{algoIndex}.fHandle(...
                            flow, 1, camera.f, Algo{algoIndex});

% *************************************************************************
% Show image flow field.
% *************************************************************************
if showFlow,
    figure; 
    quiver(flow.X*10^3,flow.Y*10^3,flow.Dx*10^3,flow.Dy*10^3,'-k'); 
    axis([-6 6 -6 6]);
    axis square; 
    xlabel('x (mm)');
    ylabel('y (mm)');
end

% *************************************************************************
% Print the grond-truth and estimated self-motion in the console.
% *************************************************************************
fprintf('Ground-truth self-motion\n');
motion.fPrint(motion);
fprintf('\nEstimated self-motion\n');
motionEst.fPrint(motionEst);

logger
