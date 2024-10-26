% This file is part of OctCLA, Copyright (c) 2022 Josh Fogg, released
% under the MIT License. See: https://github.com/Foggalong/OctCLA

% This particular example gives an indication of CLA's sensitivity to
% inputs. For a given problem, we examine how the Pareto frontier's
% turning points change when faced with small differences in expected
% return and (separately) small differences in covariance.

% set output formatting
format compact
format long

% add implementation functions to the path
addpath(genpath("../finance/"))
addpath(genpath("../genetics/"))

% define problem variables that will be consistent throughout
S  = [1];  % MATLAB will make a double, but that's fine
D  = [2,3];
lb = [0.0; 0.0; 0.0];
ub = [1.0; 1.0; 1.0];

% first keep the covariance constant, vary expected return
covar = diag([1, 5, 10]);
mu1 = [1; 1.1; 1];
sols1 = calculate_turningpoints_gen(mu1, covar, lb, ub, S, D)'
mu2 = [1; 1; 1.1];
sols2 = calculate_turningpoints_gen(mu2, covar, lb, ub, S, D)'

% now keep expected return constant, but vary covariance
mu = [1; 5; 3];
covar3 = diag([1, 5, 5.1]);
sols3 = calculate_turningpoints_gen(mu, covar3, lb, ub, S, D)'
covar4 = diag([1, 5.1, 5]);
sols4 = calculate_turningpoints_gen(mu, covar4, lb, ub, S, D)'
