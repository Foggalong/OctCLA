% This file is part of OctCLA, Copyright (c) 2022 Josh Fogg, released
% under the MIT License. See: https://github.com/Foggalong/OctCLA

% Below we demonstrate using OctCLA to solve a traditional portfolio
% optimization problem by loading data from a file, followed by an
% optimal contribution selection problem where input data is defined
% through in-code variables. The computed turning points are checked
% both using the KKT conditions and against the known solutions. 

% set output formatting
format compact
format long


disp("FINANCE CODE DEMO")
addpath(genpath("finance/"))  % need functions in path

% read data from input file into variables
data = csvread("demos/data/input.csv");
mu = data(2,:)';
lb = data(3,:)';
ub = data(4,:)';
covar = data(5:end,:);

% calculate turning points using OctCLA finance functions
tic
sols = calculate_turningpoints(mu, covar, lb, ub, 3)';
toc

% write output to file for reference
csvwrite("demos/data/output.csv", sols)

% read truth into matrix
truth = csvread("demos/data/truth.csv");

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
addpath(genpath("genetics/"))  % need functions in path

% set up the problem manually using variables
S  = [1,2];
D  = [3,4];
mu = [4; 1; 3; 2];
lb = repmat(0.0, 4, 1);
ub = repmat(1.0, 4, 1);
covar = [
    2, 0, 0, 0;
    0, 1, 0, 0;
    0, 0, 2, 0;
    0, 0, 0, 1;
];

% true turning points of this problem
truth = [
   0.5, 0.0, 0.5, 0.0;
   0.5, 0.0, 5/18, 2/9;
   1/6, 1/3, 1/3, 1/6
]';

% calculate turning points using OctCLA genetics functions
tic
sols = calculate_turningpoints_gen(mu, covar, lb, ub, S, D, 3)';
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
