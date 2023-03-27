% This file is part of OctCLA, Copyright (c) 2022 Josh Fogg, released
% under the MIT License. See: https://github.com/Foggalong/OctCLA

% TODO added a B argument which isn't in the original, need to check if it can be removed
function [ins, lam_ins, gam_ins, del_ins, b_ins, d] = move_to_bound_gen(mu, covar, invcovarF, lb, ub, F, B, S, D, lam_current, w, KKT)
   % MOVE_TO_BOUND_GEN handle the genetics CLA case where an asset moves to its bound
    %
    % As input, takes a vector (mu) of expected returns, a covariance matrix (covar),
    % the inverse of that covar when restricted to free assets (invcovarF), a vector
    % of lower bounds on assest weights (lb), a vector of upper bounds on asset
    % weights (ub), index sets for the free (F) and bounded (B) asset, index sets for
    % the sires (S) and dams (D) the current, lambda (lam_current), the current weight
    % vector (w), and the type of KKT check being carried out (KKT).
    %
    % As output, returns the i (ins), lambda (lam_ins), gamma (gam_ins), delta (del_ins),
    % and bound (b_ins) of the asset which would move to its bound. If running a full
    % KKT check, also returns the grad vector (d) which indicates which bounds each
    % asset is moving towards.
    %
    % See also, CALCULATE_TURNINGPOINTS_GEN

    % a sole free sire and sole free dam cannot move to bound
    F_D = subindex(D, F);
    F_S = subindex(S, F);
    if (length(F_D) == 1) && (length(F_S) == 1)
        ins = NaN; b_ins = NaN; gam_ins = NaN; del_ins = NaN; d = NaN; lam_ins = -inf;
        return
    end

    % pre-allocate storage vectors
    b   = zeros(length(mu), 1);  % bounds being moved towards
    lam = zeros(length(mu), 1);  % holds potential lambda values
    gam = zeros(length(mu), 1);  % gamma vector
    del = zeros(length(mu), 1);  % delta vector

    % BUG resolve duplicated variable name (b, bound vector, b linear system value)
    % probably easier at this point to rename the linear system variables since I've
    % only used them in a limited set of situations so far.

    muF          = mu(F);
    invcovarFmuF = invcovarF*muF;
    covarFB      = covar(F,B);

    % calculate derivative and multiplier updates using function
    % [gam_ins, del_ins, C, lam_new] = multiplier_update(F, B, S, D, w, invcovarF, covarFB, muF, lam_current, lb, ub, true);
    % trying to tread gamma and delta from here as vectors
    [gam_vec, del_vec, C, lam_vec] = multiplier_update(F, B, S, D, w, invcovarF, covarFB, muF, lam_current, lb, ub, true);

    % TODO check if there's a more efficient way to do this allocation
    % QUESTION is this actually even needed anymore? 
    for j = 1:length(F)
        i = F(j);
        lam(i) = lam_vec(j);
        del(i) = del_vec(j);
        gam(i) = gam_vec(j);
    end  

    % BUG need to pass b back from update, not passed currently
    % only need d if running the full KKT check
    if (KKT == 1) || (KKT == 3)
        d = b;
        d(B) = w(B);
    else
        d = NaN;
    end
    
    % check whether found new turning point
    [lam_ins, ins] = max_bounded(lam, lam_current);

    if isnan(ins)
        % other variables set to NaN/inf by max_bounded
        b_ins = NaN; gam_ins = NaN; del_ins = NaN;
    elseif (length(F_D) == 1) && (ismember(ins, D))
        % can't move sole free dam to bound
        ins = NaN; b_ins = NaN; gam_ins = NaN; del_ins = NaN; d = NaN; lam_ins = -inf;
    elseif (length(F_S) == 1) && (ismember(ins, S))
        % can't move sole free sire to bound
        ins = NaN; b_ins = NaN; gam_ins = NaN; del_ins = NaN; d = NaN; lam_ins = -inf;
    else
        % found a turning point, return b
        b_ins = b(ins);
        gam_ins = gam(ins);
        del_ins = del(ins);
    end
end
