% This file is part of OctCLA, Copyright (c) 2022 Josh Fogg, released
% under the MIT License. See: https://github.com/Foggalong/OctCLA

function result = kkt_partitioned(covar, w, lam, gam, mu, B, F, tol)
    % KKT_PARTITIONED equivalent KKT checks for constrained CLA
    %
    % Takes a lambda (lam), a gamma (gam), and (F)ree and (B)ounded
    % weight lists associated with a solution (w) to an unconstrained
    % portfolio optimisation problem, which is the equivalent problem
    % to the constrained problem of KKT_FULL, and then returns true if
    % that solution passes the KKT conditions and false otherwise.
    % Takes an optional float input which is the tollerance to which
    % the conditions are verified.
    %
    % See also, KKT_FULL

    % set default value for tolerance
    if (nargin < 8); tol = 1e-10; end  % TODO check this is sensible

    % Check KKT(2.2): weights sum to one
    if sum(w) - 1 > tol
        disp('Broke KKT(2.2)');
        w
        result = false; return
    end

    % Check KKT(1): Lagrangian is zero
    onesF = ones(size(mu(F)));
    delwL = covar(F,F)*w(F) + covar(F,B)*w(B) - lam*mu(F) - gam*onesF;
    if abs(delwL) > tol
        disp('Broke KKT(1)');
        delwL
        result = false; return
    end

    result = true;
end
