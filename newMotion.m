function motion = newMotion(Vel, Omega)
% newMotion
%   Vel     - Linear velocity in 3D space defined by three components.
%             Units are meter/frame.
%   Omega   - Rotational velocity in 3D space (ptich, yaw, roll).
%             Units are RAD/frame.
%
% RETURN
%   motion  - Structure with the above arguments as field.
%
%   Copyright (C) 2012  Florian Raudies, 04/12/2012, Boston University.
%   License, GNU GPL, free software, without any warranty.
%

% Set fields.
motion.Vel      = Vel;
motion.Omega    = Omega;
motion.fPrint   = @printMotion;

% Method for printing motions to the console.
function printMotion(motion)
fprintf('Vel=(%e, %e, %e)m/frame\n',...
    motion.Vel(1),motion.Vel(2),motion.Vel(3));
fprintf('Omega=(%e, %e, %e)deg/frame\n',...
    motion.Omega(1)*180/pi,motion.Omega(2)*180/pi,motion.Omega(3)*180/pi);

