% This file is part of OctCLA, Copyright (c) 2022 Josh Fogg, released
% under the MIT License. See: https://github.com/Foggalong/OctCLA

function [F, B, w] = starting_solution_gen(mu, lb, ub, S, D, tol)
    % STARTING_SOLUTION_GEN return starting solution for genetics CLA 
    %
    % Takes a vector of expected returns (mu), a vector of lower
    % bounds on assest weights (lb), a vector of upper bounds on
    % asset weights (ub), an index set for sires (S), and an index
    % set for dams (D) as input, then returns the starting solution
    % for CLA in the form of an index list of (F)ree and (B)ounded
    % assests and starting weight vector (w).
    %
    % See also, STARTING_SOLUTION, CALCULATE_TURNINGPOINTS_GEN

    % set default value for tolerance
    if (nargin < 6); tol = 1e-10; end  % TODO check this is sensible

    % start with all assets on their lower bound
    w = lb;

    % starting with inf ensures first value satisfies mi(i) < mu_max
    mu_max = inf;
    % increase sire weights in descending order of expected return
    while (abs(sum(w(S)) - 0.5) > tol)
        [mu_max, i] = max_conditional(mu, mu_max, S);
        w(i) = min(ub(i), lb(i)+0.5-sum(w(S)));
    end
    % only one sire starts free
    F = [i];
    B = S(S ~= i);

    % startting with inf ensures first value satisfies mi(i) < mu_max
    mu_max = inf;
    % increase dam weights in descending order of expected return
    while (abs(sum(w(D)) - 0.5) > tol)
        [mu_max, i] = max_conditional(mu, mu_max, D);
        w(i) = min(ub(i), lb(i)+0.5-sum(w(D)));
    end
    % only one dam starts free
    F = [F, i];
    % all other dams start on bounds
    B = [B, D(D ~= i)];
end
