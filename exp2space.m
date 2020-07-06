function X = exp2space(minVal, maxVal, num)
% exp2space
%   minVal     - Min value of interval.
%   maxVal     - Max value of interval.
%   num        - Samples in interval.
%
% RETURNS
%   X          - X = 2.^linspace(minVal,maxVal,num).
%
%   Copyright (C) 2012  Florian Raudies, 04/12/2012, Boston University.
%   License, GNU GPL, free software, without any warranty.
%

X = 2.^linspace(minVal,maxVal,num);