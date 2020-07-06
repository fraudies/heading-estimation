function noise = newGaussianNoise(sigma, mu)
% newGaussianNoise
%   sigma   - Standard deviation of Gaussian noise.
%   mu      - Mean of Gaussian noise.
%
% RETURN
%   noise   - Structure with the above parameters as arguments.
%
% DESCRIPTION
%   This noise can be applied to a flow field using the function handle of 
%   the structure 'fHandle' which computes the following expression:
%
%       Dx = Dx + N(sigma,mu)
%       Dy = Dy + N(sigma,mu)
%
%   with N being a normal distribution function of standard deviation sigma
%   and mean value mu.
%
%   Copyright (C) 2012  Florian Raudies, 04/12/2012, Boston University.
%   License, GNU GPL, free software, without any warranty.

noise.sigma     = sigma;
noise.mu        = mu;
noise.fHandle   = @addGaussianNoise2Flow;


% *************************************************************************
% Private function that adds Gaussian noise to the image flow.
% *************************************************************************
function flow = addGaussianNoise2Flow(noise, flow)
% addGaussianNoise2Flow
%   noise   - Structure that describes the noise. Mandatory fields are:
%             'sigma' and 'mu'.
%   flow    - Structure that contains the flow. Mandatory fields are: 'Dx'
%             and 'Dy'.
%
% RETURN
%   flow    - Same structure input with 'Dx' and 'Dy' manipulated.
%
Dim = size(flow.Dx);
flow.Dx = flow.Dx + noise.mu + noise.sigma*randn(Dim);
flow.Dy = flow.Dy + noise.mu + noise.sigma*randn(Dim);
