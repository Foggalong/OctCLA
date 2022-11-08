% This file is part of OctCLA, Copyright (c) 2022 Josh Fogg, released
% under the MIT License. See: https://github.com/Foggalong/OctCLA

% set output formatting
format compact
format long

disp("FINANCE CODE DEMO")
% TODO write a suitable description of this code
addpath(genpath("finance/"))

% read data from input file into variables
data = csvread("demos/input.csv");
mu = data(2,:)';
lb = data(3,:)';
ub = data(4,:)';
covar = data(5:end,:);

% calculate turning points using OctCLA functions
tic
sols = calculate_turningpoints(mu, covar, lb, ub, 3)';
toc

% write output to file for reference
csvwrite("demos/output.csv", sols)

% read truth into matrix
truth = csvread("demos/truth.csv");

% check if any entry in absolute error matrix greater than tolerance
tol = 1e-10;
abs_error = abs(truth - sols);
if any(any(abs_error > tol) > 0) then
    disp("ERROR!")
    abs_error
else
    disp("Success!")
end

disp("")  % padding

disp("GENETICS CODE DEMO")
% TODO write a suitable description of this code
addpath(genpath("genetics/"))

% set up the problem manually using variables
S  = [1];  % MATLAB will make a double, but that's fine
D  = [2,3];
mu = [1; 5; 2];
lb = [0.0; 0.0; 0.0];
ub = [1.0; 1.0; 1.0];
covar = [
    1, 0, 0;
    0, 5, 0;
    0, 0, 3;
];

% true turning points of this problem
truth = [
    0.50000, 0.50000, 0.00000;
    0.50000, 0.18750, 0.31250
]';

% calculate turning points using OctCLA genetics functions
tic
sols = calculate_turningpoints_gen(mu, covar, lb, ub, S, D)';
toc

% return turning points to terminal
sols

% check if any entry in absolute error matrix greater than tolerance
tol = 1e-10;
abs_error = abs(truth - sols);
if any(any(abs_error > tol) > 0) then
    disp("ERROR!")
    abs_error
else
    disp("Success!")
end

