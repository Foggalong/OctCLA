% This file is part of OctCLA, Copyright (c) 2022 Josh Fogg, released
% under the MIT License. See: https://github.com/Foggalong/OctCLA

function [F, B, w] = starting_solution(mu, lb, ub)
    % STARTING_SOLUTION return starting solution for CLA
    %
    % Takes a vector of expected returns (mu), a vector of lower
    % bounds on assest weights (lb), and a vector of upper bounds on
    % asset weights (ub) as input, then returns the starting solution
    % for CLA in the form of an index list of (F)ree and (B)ounded
    % assests and starting weight vector (w).
    %
    % See also, CALCULATE_TURNINGPOINTS

    % start with all assets on their lower bound
    w = lb;
    % set to inf to make first max_bounded equivalent to max(mu)
    mu_max = inf;
    % increase assest weights in descending order of expected return
    while sum(w) < 1
        % move to the next highest mu(i) after mu_max
        [mu_max, i] = max_bounded(mu, mu_max);
        % set to maximum or use remaining capacity, whichever's lower
        w(i) = min(ub(i), lb(i)+1-sum(w));
    end
    % only one asset starts free
    F = [i];
    % all other assets start on bounds
    B = 1:length(mu);
    B = B(B ~= i);
end
