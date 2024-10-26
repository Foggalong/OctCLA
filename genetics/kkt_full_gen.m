% This file is part of OctCLA, Copyright (c) 2022 Josh Fogg, released
% under the MIT License. See: https://github.com/Foggalong/OctCLA

% TODO implement replacement function for KKT_FULL. Did this once but
% seem to have lost it at some point :(

function result = kkt_full_gen(lb, ub, d, w_t, covar, lam, mu, gam, del, S, D, tol)
    % TODO write a docstring

    % set default value for tolerance
    if (nargin < 12); tol = 1e-10; end  % TODO check this is sensible

    % Check KKT(2.2.s): sire weights sum to half
    if sum(w_t(S)) - 0.5 > tol
        disp('Broke KKT(2.2.s)');
        w_t_S = w_t(S)
        result = false; return
    end

    % Check KKT(2.2.d): dam weights sum to half
    if sum(w_t(D)) - 0.5 > tol
        disp('Broke KKT(2.2.d)');
        w_t_D = w_t(D)
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

    % vector(i) has gamma if i Sire, delta if i Dam  
    gam_del_vector = zeros(size(mu));
    gam_del_vector(S) = gam;
    gam_del_vector(D) = del;

    % KKT(1) is assumed through this statement
    tau = C*(covar*w_t - lam*mu - gam_del_vector);

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
