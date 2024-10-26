% This file is part of OctCLA, Copyright (c) 2022 Josh Fogg, released
% under the MIT License. See: https://github.com/Foggalong/OctCLA

function [outs, lam_outs, d] = becomes_free(mu, covar, invcovarF, lb, ub, F, B, lam_current, w, KKT)
    % BECOMES_FREE handle the CLA case where an asset becomes free

    % As input, takes a vector (mu) of expected returns, a covariance
    % matrix (covar), the inverse of that covar when restricted to
    % free assets (invcovarF), a vector of lower bounds on assest
    % weights (lb), a vector of upper bounds on asset weights (ub),
    % index sets for the free (F) and bounded (B) asset, the current
    % lambda (lam_current), the current weight vector (w), and the
    % type of KKT check being carried out (KKT).
    %
    % As output, returns the i (outs) and lambda (lam_outs) of the
    % asset which moves away from its bound (i.e. which becomes free).
    % If running a full KKT check, also returns the grad vector (d)
    % which indicates which bounds each asset is moving towards.
    %
    % See also, CALCULATE_TURNINGPOINTS

    % TODO added a B argument not in original, check if it can be removed

    % skip proceedure if all assets are free
    if (length(F) == length(mu))
        outs = NaN; d = NaN; lam_outs = -inf;
        return
    end

    lam = zeros(length(mu), 1);  % lam_i = lambda if i becomes free
    % only need D if running the full KKT check
    if (KKT == 1) || (KKT == 3)
        % matrix of potential d vectors
        possible_d = zeros(length(mu), length(mu));
    end

    for i = B
        % update the inverse to reflect additional row and column
        invcovarFi = inverse_grow(invcovarF, covar(F,i), covar(i,i));
        % free weight i
        Fi = [F, i];     % F = Fu{i}
        Bi = B(B ~= i);  % B = B\{i}
        % need to index outer matrix, as per NOTE A1
        j = length(Fi);  % Fi[j] = i; i last element in Fi by above

        % calculate derivative  % TODO tidy this up
        Ci1 = -sum(sum(invcovarFi))*(invcovarFi*mu(Fi));
        Ci2 = (sum(invcovarFi)*mu(Fi))*sum(invcovarFi, 2);
        Ci = Ci1 + Ci2;

        % if running full KKT check, save d vector
        if (KKT == 1) || (KKT == 3)
            for l = 1:length(Fi)
                k = Fi(l);
                if Ci(l) > 0; possible_d(k, i) = ub(l); end
                if Ci(l) < 0; possible_d(k, i) = lb(l); end
            end
            possible_d(Bi, i) = w(Bi);
        end

        % handle case in NOTE A2
        lami_p1 = sum(invcovarFi(j,:), 2);
        if isempty(Bi)
            lami_p2 = sum(sum(invcovarFi))*w(i);
        else
            % calculate lambda using shortcuts  % TODO tidy this up
            lami_p1 = (1-sum(w(Bi))+sum(invcovarFi)*(covar(Fi,Bi)*w(Bi))) * lami_p1;
            lami_p2_q2 = invcovarFi*(covar(Fi,Bi)*w(Bi));
            lami_p2 = sum(sum(invcovarFi))*(w(i)+lami_p2_q2(j));
        end

        lam(i) = (lami_p1-lami_p2)/Ci(j);
    end

    % check whether found new turning point
    [lam_outs, outs] = max_bounded(lam, lam_current);

    % only have d if outs defined and doing full KKT check
    if isnan(outs) || (KKT == 0) || (KKT == 2)  
        d = NaN;
    else
        d = possible_d(:, outs);
    end
end
