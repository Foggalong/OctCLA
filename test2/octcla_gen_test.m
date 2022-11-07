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

% This has turning points
% ws =
%    0.50000   0.50000   0.00000
%    0.50000   0.18750   0.31250

% source('../octcla.m')      % only Octave compatible
% source('../octcla-gen.m')  % only Octave compatible
tic
sols = calculate_turningpoints_gen(mu, covar, lb, ub, S, D)'
toc
