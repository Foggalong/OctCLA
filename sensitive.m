% This file is part of OctCLA, Copyright (c) 2022 Josh Fogg, released
% under the MIT License. See: https://github.com/Foggalong/OctCLA

% set output formatting
format compact
format long

% TODO write a suitable description of this code
addpath(genpath("finance/"))
addpath(genpath("genetics/"))

% set up the problem manually using variables
S  = [1];  % MATLAB will make a double, but that's fine
D  = [2,3];
mu = [1; 1.1; 1];
lb = [0.0; 0.0; 0.0];
ub = [1.0; 1.0; 1.0];
covar = [
    1, 0, 0;
    0, 5, 0;
    0, 0, 5;
];

% calculate turning points using OctCLA genetics functions
tic
sols = calculate_turningpoints_gen(mu, covar, lb, ub, S, D)';
toc

% return turning points to terminal
sols

mu = [1; 1; 1.1];
covar = [
    1, 0, 0;
    0, 5, 0;
    0, 0, 5;
];

% calculate turning points using OctCLA genetics functions
tic
sols = calculate_turningpoints_gen(mu, covar, lb, ub, S, D)';
toc

% return turning points to terminal
sols
