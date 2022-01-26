% This octcla-test.m is part of OctCLA, Copyright (c) 2022 Josh Fogg, released
% under the MIT License. For more info see https://github.com/Foggalong/OctCLA

% set output formatting
format compact
format long

% read data from input file into variables
data = csvread('input.csv');
mu = data(2,:)';
lb = data(3,:)';
ub = data(4,:)';
covar = data(5:end,:);

% calculate turning points using OctCLA functions
source('../octcla.m')
sols = calculate_turningpoints(mu, covar, lb, ub, KKT=3)';

% write output to file for reference
csvwrite('output.csv', sols)

% read truth into matrix
truth = csvread('truth.csv');

% check if any entry in absolute error matrix greater than tolerance
tol = 1e-10;
abs_error = abs(truth - sols);
if any(any(abs_error > tol) > 0) then
    printf("ERROR!\n")
    abs_error
else
    printf("Success!\n")
end
