function flow = newImgFlow(X, Y, Dx, Dy)
% newImgFlow
%   X       - x-coordinate in the image plane.
%   Y       - y-coordinate.
%   Dx      - x-component of the flow vector.
%   Dy      - y-component.
%
% RETURN
%   flow    - Structure with the above arguments as fields.
%
%   Copyright (C) 2012  Florian Raudies, 04/12/2012, Boston University.
%   License, GNU GPL, free software, without any warranty.
%

flow.X  = X;
flow.Y  = Y;
flow.Dx = Dx;
flow.Dy = Dy;