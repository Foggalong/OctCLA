% This file is part of OctCLA, Copyright (c) 2022 Josh Fogg, released
% under the MIT License. See: https://github.com/Foggalong/OctCLA

function [x_max, i_max] = max_bounded(x, x_bound, tol)
    % MAX_BOUNDED return highest value and index less than a bound
    %
    % Takes a vector x and a number x_bound as inputs and then returns
    % the value x_max and index i_max of the highest value in vector x 
    % satisfying the condition x(i_max) < x_bound. If x_bound = inf it
    % is equivalent to max(x). Can also take an optional value tol
    % which specifies them tolerance (default: 10^-10).
    %
    % See also, MAX.

    % set default value for tolerance
    if (nargin < 3); tol = 1e-10; end

    % starting values catch the case where x empty
    x_max = -inf;
    i_max = 0;

    % want the index, so iterate over location rather than value
    for i = 1:length(x)
        % if x_bound=inf, second condition always satisfied
        if (tol < x(i)-x_max && tol < x_bound-x(i))
            x_max = x(i); i_max = i;
        end
    end

    % if no max found, return NaN
    if i_max == 0
        x_max = NaN; i_max = NaN;
    end
end
