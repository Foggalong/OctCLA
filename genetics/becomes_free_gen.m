% This file is part of OctCLA, Copyright (c) 2022 Josh Fogg, released
% under the MIT License. See: https://github.com/Foggalong/OctCLA

% TODO added a B argument which isn't in the original, need to check if it can be removed
function [outs, lam_outs, gam_outs, del_outs d] = becomes_free_gen(mu, covar, invcovarF, muF, lb, ub, F, B, S, D, lam_current, w, KKT)
    % BECOMES_FREE_GEN handle the genetics CLA case where an asset becomes free
    %
    % As input, takes a vector (mu) of expected returns, a covariance matrix (covar),
    % the inverse of that covar when restricted to free assets (invcovarF), a vector
    % of lower bounds on assest weights (lb), a vector of upper bounds on asset
    % weights (ub), index sets for the free (F) and bounded (B) asset, index sets for
    % the sires (S) and dams (D), the current lambda (lam_current), the current weight
    % vector (w), and the type of KKT check being carried out (KKT).
    %
    % As output, returns the i (outs), lambda (lam_outs), gamma (gam_outs), and delta
    % (del_outs) of the asset which moves away from its bound (i.e. which becomes free).
    % If running a full KKT check, also returns the grad vector (d) which indicates
    % which bounds each asset is moving towards.
    %
    % See also, CALCULATE_TURNINGPOINTS_GEN

    % skip proceedure if all assets are free
    if (length(F) == length(mu))
        outs = NaN; gam_outs = NaN; del_outs = NaN; d = NaN; lam_outs = -inf;
        return
    end

    lam = zeros(length(mu), 1);  % lambda vector
    gam = zeros(length(mu), 1);  % gamma vector
    del = zeros(length(mu), 1);  % delta vector

    % only need D if running the full KKT check
    if (KKT == 1) || (KKT == 3)
        possible_d = zeros(length(mu), length(mu));  % matrix of potential d vectors
    end

    for i = B
        % update the inverse
        a = covar(F, i);  % BUG another duplicated variable name
        alpha = covar(i, i);
        invcovarFi = inverse_grow(invcovarF, a, alpha);
        % free weight i
        Fi = [F, i];     % F = Fu{i}
        Bi = B(B ~= i);  % B = B\{i}
        % don't need to update S/D, constant between iterations

        % additional shortcuts needed just in genetics version
        covarFiBi = covar(Fi,Bi);

        % need to index outer matrix, as per NOTE A1
        j = length(Fi);  % Fi[j] = i; i last element in Fi by construction

        % calculate derivative and multiplier updates using function
        % BUG this is calculating MANY unneeded values just to get the ith
        [gam_vec, del_vec, Ci, lam_vec] = multiplier_update(Fi, Bi, S, D, w, invcovarFi, covarFiBi, mu(Fi), lam_current, lb, ub, false);
        lam(i) = lam_vec(j);
        del(i) = del_vec(j);
        gam(i) = gam_vec(j);
        % if running full KKT check, save d vector
        if (KKT == 1) || (KKT == 3)
            for l = 1:length(Fi)
                k = Fi(l);
                if Ci(l) > 0; possible_d(k, i) = ub(l); end
                if Ci(l) < 0; possible_d(k, i) = lb(l); end
            end
            possible_d(Bi, i) = w(Bi);  % TODO factor this into a KKT if statement
        end
    end

    % check whether found new turning point
    [lam_outs, outs] = max_bounded(lam, lam_current);

    % need to check rather than assume or will get an index error
    if (outs ~= NaN)
        % select correponding gamma and delta multipliers
        gam_outs = gam(outs);
        del_outs = del(outs);
    else
        gam_outs = NaN; del_outs = NaN;
    end

    % only have d if outs defined and doing full KKT check
    if (outs == NaN) || (KKT == 0) || (KKT == 2)  
        d = NaN;
    else
        d = possible_d(:, outs);
    end
end
