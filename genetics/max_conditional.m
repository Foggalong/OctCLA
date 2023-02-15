% This file is part of OctCLA, Copyright (c) 2022 Josh Fogg, released
% under the MIT License. See: https://github.com/Foggalong/OctCLA

function [x_max, i_max] = max_conditional(x, x_bound, index_set, tol)
    % MAX_CONDITIONAL return value and index of max item meeting conds
    %
    % Takes a vector x, a number x_bound, and a set of indices
    % index_set as inputs and returns the value x_max and index i_max
    % of the highest value in x satisfying conditions that i_max is in
    % index_set and x(i_max) < x_bound. Note if index_set is
    % 1:length(x) the function is equivalent to max_bounded. Takes an
    % optional value tol which specifies tolerance (default: 10^-10).
    %
    % See also, MAX, MAX_BOUNDED

    % set default value for tolerance
    if (nargin < 4); tol = 1e-10; end

    % starting values catch the case where x empty
    x_max = -inf;
    i_max = 0;

    % only interested in indicies from index_set so iterate over them
    for i = index_set
        % if x_bound=inf second check always satified
        if (tol < x(i)-x_max && tol < x_bound-x(i))
            x_max = x(i); i_max = i;
        end
    end

    % if no max found, return NaN
    if i_max == 0
        x_max = NaN; i_max = NaN;
    end
end
