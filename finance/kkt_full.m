% This file is part of OctCLA, Copyright (c) 2022 Josh Fogg, released
% under the MIT License. See: https://github.com/Foggalong/OctCLA

function result = kkt_full(lb, ub, d, w_t, covar, lam, mu, gam, tol)
    % KKT_FULL perform a full KKT check on a constrained CLA problem
    %
    % Takes lambda (lam), gamma (gam), and derivative vector (d)
    % values associated with a solution (w_t) to a portfolio
    % optimisation problem with given covariance matrix (covar),
    % expected return vector (mu), and bound constraints (lb, ub) and
    % then returns true if that solution passes the KKT conditions and
    % false otherwise. Takes an optional float input which is the
    % tollerance to which the conditions are verified.
    %
    % See also, KKT_PARTITIONED

    % set default value for tolerance
    if (nargin < 9); tol = 1e-10; end  % TODO check this is sensible

    % Check KKT(2.2): weights sum to one
    if sum(w_t) - 1 > tol
        disp('Broke KKT(2.2)');
        w_t = w_t
        result = false; return
    end

    % The ith entry of vector d contains the value of the bound which
    % the ith weight was moving toward in the solution. Thus the
    % diagonal matrix C is such that the (i,i)th entry is +1 if the
    % ith weight was moving toward the lower bound and -1 if it was
    % moving toward the upper bound.
    C = diag([d==lb] - [d==ub]);

    % Check KKT(3): weights within bounds
    if C*w_t-d > tol
        disp('Broke KKT(3)'); 
        d_Cw = d-C*w_t
        result = false; return
    end

    % KKT(1) is assumed through this statement
    tau = C*(covar*w_t - lam*mu - gam*ones(size(mu)));

    % Check KKT(4): Tau positive
    if tau < -tol
        disp('Broke KKT(4)');
        tau = tau
        result = false; return
    end

    % Check KKT(5): tau'*(Cw-d) = 0
    if abs(tau'*(C*w_t - d)) > tol
        disp('Broke KKT(5)');
        d = d
        Cw_d = C*w_t-d
        tau = tau
        tCw_d = tau'*(C*w_t - d)
        result = false; return
    end

    result = true;
end
