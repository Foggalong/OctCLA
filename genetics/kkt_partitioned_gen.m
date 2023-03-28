% This file is part of OctCLA, Copyright (c) 2022 Josh Fogg, released
% under the MIT License. See: https://github.com/Foggalong/OctCLA

function result = kkt_partitioned_gen(covar, w, lam, gam, del, mu, B, F, S, D, tol)
    % KKT_PARTITIONED_GEN equivalent  $ TODO update
    %
    % TODO write a docstring
    %
    % See also, KKT_FULL_GEN, KKT_PARTITIONED

    % set default value for tolerance
    if (nargin < 11); tol = 1e-10; end  % TODO check this is sensible

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

    % vector(i) has gamma if i Sire, delta if i Dam  
    gam_del_vector = zeros(size(mu));
    gam_del_vector(S) = gam;
    gam_del_vector(D) = del;

    % Check KKT(1): Lagrangian is zero
    delwL = covar(F,F)*w(F) + covar(F,B)*w(B) - lam*mu(F) - gam_del_vector;
    if abs(delwL) > tol
        disp('Broke KKT(1)');
        delwL
        result = false; return
    end

    result = true;
end
