% This file is part of OctCLA, Copyright (c) 2022 Josh Fogg, released
% under the MIT License. See: https://github.com/Foggalong/OctCLA

function [ins, lam_ins, b_ins, d] = move_to_bound(mu, covar, invcovarF, lb, ub, F, B, lam_current, w, KKT)
    % MOVE_TO_BOUND handle CLA case where an asset moves to its bound
    %
    % As input, takes a vector (mu) of expected returns, a covariance
    % matrix (covar), the inverse of that covar when restricted to
    % free assets (invcovarF), a vector of lower bounds on assest
    % weights (lb), a vector of upper bounds on asset weights (ub),
    % index sets for the free (F) and bounded (B) asset, the current
    % lambda (lam_current), the current weight vector (w), and the
    % type of KKT check being carried out (KKT).
    %
    % As output, returns the i (ins), lambda (lam_ins), and bound
    % (b_ins) of the asset which would move to its bound. If running a
    % full KKT check, also returns the grad vector (d) which indicates
    % which bounds each asset is moving towards.
    %
    % See also, CALCULATE_TURNINGPOINTS

    % TODO added a B argument not in original, check if it can be removed

    % a sole free asset cannot move to bound
    if length(F) == 1
        ins = NaN; b_ins = NaN; d = NaN; lam_ins = -inf;
        return
    end

    b   = zeros(length(mu), 1);  % b vector
    lam = zeros(length(mu), 1);  % lambda vector

    muF          = mu(F);
    invcovarFmuF = invcovarF*muF;
    covarFB      = covar(F,B);

    % calculate derivative, using shortcuts
    C = -sum(sum(invcovarF))*(invcovarFmuF) + sum(invcovarFmuF)*sum(invcovarF,2);

    % precalculate lam_p1 if it will be needed
    if ~isempty(B)
        lam_p1 = (1-sum(w(B))+sum(invcovarF)*(covarFB*w(B)))*sum(invcovarF,2);
    else
        lam_p1 = sum(invcovarF, 2);
    end

    for j = 1:length(F)
        % need to index outer matrix, as per NOTE A1
        i = F(j);
        if C(j) > 0
            b(i) = ub(i);
        elseif C(j) < 0
            b(i) = lb(i);
        else  % C(j) == 0
            % since |F| > 1, C(j) == 0 iff all mu(i) equal
            continue  % TODO check how properly to handle this
        end

        % calculate lambda using shortcuts  % TODO tidy this up some
        if isempty(B)
            lami_p2 = sum(sum(invcovarF))*b(i);
        else
            invcovarFcovarFBwB = invcovarF*(covarFB*w(B));
            lami_p2 = sum(sum(invcovarF))*(b(i)+invcovarFcovarFBwB(j));
        end

        lam(i) = (lam_p1(j)-lami_p2)/C(j);
    end

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
        b_ins = NaN;
    else
        b_ins = b(ins);
    end
end