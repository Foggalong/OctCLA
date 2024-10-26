% This file is part of OctCLA, Copyright (c) 2022 Josh Fogg, released
% under the MIT License. See: https://github.com/Foggalong/OctCLA

function [gam, del, C, lambda, bound] = multiplier_update(F, B, S, D, w, invcovarF, covarFB, muF, lam_current, lb, ub, m2b, tol)
    % MULTIPLIER_UPDATE handle the updating of Lagrangian multipliers
    %
    % To avoid passing multiple variables between other functions, given the
    % more complicated update equations in the genetics context it turns out
    % to update our multipliers through a single function. This avoids having
    % lots of duplicate code split across multiple places, but comes at the
    % cost of being longer than ideal.

    % As input, it takes the index sets for the free (F) and bounded (B) assets,
    % the index sets of sires (S) and dams (D), the current portfolio (w), the
    % current `invcovarF`, `covarFB`, and `muF` (avoids duplicate calculation),
    % the current value of lambda (lam_current), lower bound on weights (lb),
    % and upper bound on weights (ub). The final mandatory argument `m2b' 
    % determines which branch of the algorithm the function is being called in.
    % If m2b = true, it's in move_to_bound, otherwise it's in becomes_free. This
    % changes how the `bound` vector is calculated.
    %
    % As output, returns the updated gamma (gam), delta (del), derivative (C),
    % `lambda`, and the bound each variable is moving towards (bound). These
    % correspond to the updated multipliers for the given inputs.
    %
    % See also, MOVE_TO_BOUND_GEN, BECOMES_FREE_GEN, CALCULATE_TURNINGPOINTS_GEN

    % set default value for tolerance
    if (nargin < 13); tol = 1e-10; end  % TODO check this is sensible

    % want D and S to be indexes in free assets, not in all the assets
    F_D = subindex(D, F);
    F_S = subindex(S, F);

    % column sums of the four quadrants of invcovarF
    eF_DinvcovarF_DD = sum(invcovarF(F_D,F_D), 1);
    eF_DinvcovarF_DS = sum(invcovarF(F_D,F_S), 1);
    eF_SinvcovarF_SD = sum(invcovarF(F_S,F_D), 1);
    eF_SinvcovarF_SS = sum(invcovarF(F_S,F_S), 1);

    % matrix entries on the LHS of the multipler linear system
    a = sum(eF_DinvcovarF_DS);
    b = sum(eF_DinvcovarF_DD);
    c = sum(eF_SinvcovarF_SS);
    % d = a, so doesn't need repeating

    % determinant of the multiplier linear system
    determinant = a^2 - b*c;
    % BUG work out when this occurs; could be impossible in well formed problems
    if (abs(determinant) < tol)
        % TODO work out how to handle situations
        disp("ERROR! Linear system has no solutions")
    end


    % TODO refactor to involve code duplication
    % components of x & y are reused elsewhere so split up
    x_p1 = 0.5;
    y_p1 = 0.5;
    % if B empty we get cancellations, so only consider non-empty case
    if ~isempty(B)
        % for x & y also need sub-indexed wB
        B_D = subindex(D, B);
        B_S = subindex(S, B);
        wB = w(B);
        wB_S = wB(B_S);
        wB_D = wB(B_D);
        % first part of x and y when B non-empty
        if isempty(B_D)
            % cancellations if no bounded dams
            x_p1 = x_p1 + eF_DinvcovarF_DD*(covarFB(F_D, B_S)*wB_S);
            y_p1 = y_p1 - sum(wB_S, 1) + eF_SinvcovarF_SS*(covarFB(F_S, B_S)*wB_S);
        elseif isempty(B_S)
            % cancellations if no bounded sires
            x_p1 = x_p1 - sum(wB_D, 1) + eF_DinvcovarF_DD*(covarFB(F_D, B_D)*wB_D);
            y_p1 = y_p1 + eF_SinvcovarF_SS*(covarFB(F_S, B_D)*wB_D);
        else
            % no cancellations if we've bounded dams and sires
            x_p1 = x_p1 - sum(wB_D, 1) + eF_DinvcovarF_DD*(covarFB(F_D, B_S)*wB_S) + eF_DinvcovarF_DD*(covarFB(F_D, B_D)*wB_D);
            y_p1 = y_p1 - sum(wB_S, 1) + eF_SinvcovarF_SS*(covarFB(F_S, B_S)*wB_S) + eF_SinvcovarF_SS*(covarFB(F_S, B_D)*wB_D);
        end
    end

    x_p2 = eF_DinvcovarF_DS*muF(F_S) + eF_DinvcovarF_DD*muF(F_D);
    y_p2 = eF_SinvcovarF_SS*muF(F_S) + eF_SinvcovarF_SD*muF(F_D);

    % final multiplier update is much simpler
    if (lam_current == 0)
        gam = (a*x_p1 - b*y_p1)/determinant;
        del = (a*y_p1 - c*x_p1)/determinant;
        lambda = 0;
        % unneeded outputs in this case
        C = NaN;
        bound = NaN;
        return
    end

    % derivatives of previous multipliers with respect to lambda
    dGam_dLam = (a*x_p2 - b*y_p2)/determinant;
    dDel_dLam = (a*y_p2 - c*x_p2)/determinant;

    % C(i) indicates if asset i is moving to its upper or lower bound
    gam_del_derivatives = zeros(size(F))';
    gam_del_derivatives(F_S) = dGam_dLam;
    gam_del_derivatives(F_D) = dDel_dLam;
    C = invcovarF * (gam_del_derivatives + muF);

    % lambda update equation is separated into chunks for readability.
    % the first of these chunks is determined by whether we're considering
    % assets becoming free or moving to a bound.

    if m2b
        % <compute vector of which bound each asset is moving towards>
        bound = zeros(length(w), 1);

        % in move_to_bound, so b(i) determined based on derivative
        for j = 1:length(F)
            i = F(j);  % need to index outer matrix, as per NOTE A1
            if C(j) < 0
                bound(i) = ub(i);
            elseif C(j) > 0
                bound(i) = lb(i);
            else  % C(j) == 0
                continue % since |F| > 1, C(j) == 0 iff all mu(i) equal
            end
        end
        % precalculate lam_num_p1
        if isempty(B)
            lam_num_p1 = determinant*bound(F);
        else
            lam_num_p1 = determinant*(bound(F) + invcovarF*covarFB*w(B));
        end
    else
        bound = NaN;
        % in becomes_free, so b(i) = w(i)
        if isempty(B)
            lam_num_p1 = determinant*w(F);
        else
            lam_num_p1 = determinant*(w(F) + invcovarF*covarFB*w(B));
        end
    end

    % we might think to use `repmat` here ease, but any code like
    %     [repmat(..., size(F_S)), repmat(..., size(F_D))]'
    % will mean these lambda entry calculations become ordered into
    % sires and dams, leading subsequent calculations to be wrong.
    
    lam_num_p2 = zeros(size(F))';
    lam_num_p2(F_S) = b*y_p1 - a*x_p1;
    lam_num_p2(F_D) = c*x_p1 - a*y_p1;
    
    lam_den_p1 = zeros(size(F))';
    lam_den_p1(F_S) = b*y_p2 - a*x_p2;
    lam_den_p1(F_D) = c*x_p2 - a*y_p2;

    lam_den_p2 = invcovarF*(lam_den_p1 + determinant*muF);
    lambda = (lam_num_p1 + invcovarF*lam_num_p2)./lam_den_p2;

    % matrix entries on the RHS of the multipler linear system
    x = x_p1 - lambda*x_p2;
    y = y_p1 - lambda*y_p2;

    % update multipliers by direct calculateion of system solutions
    gam = (a*x - b*y)/determinant;
    del = (a*y - c*x)/determinant;

    % NOTE addityional x, y, gam, and del calculations are being done here
    % than necessary (because lambda is a vector and we only choose one).
    % Would be worth checking if we can make this better
end
