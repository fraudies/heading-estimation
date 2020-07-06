function logger = newEstSelfMotionLogger()
% newEstSelfMotionLogger
%   no arguments
%
% RETURNS
%   logger  - Structure with pre-set fields.
%
%   Copyright (C) 2012  Florian Raudies, 04/12/2012, Boston University.
%   License, GNU GPL, free software, without any warranty.

% For fix-point or Gauss-Newton (GN) iteration.
logger.initMean      = 0;
logger.initStd       = 0;
logger.iterMean      = 0;
logger.iterStd       = 0;
% For iterations used for m-functions.
logger.mFuncIter     = 0;