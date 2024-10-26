% This file is part of OctCLA, Copyright (c) 2022 Josh Fogg, released
% under the MIT License. See: https://github.com/Foggalong/OctCLA

function [gam, del, C, lambda] = multiplier_update(F, B, S, D, w, invcovarF, covarFB, muF, lam_current, lb, ub, m2b, tol)
    % MULTIPLIER_UPDATE handle the updating of lagrangian multipliers
    % 
    % TODO Write a proper function description.
    % TODO find a way to make ub and lb non-essential
    %
    % NOTE This started mostly me just working through the full set of multiplier
    % update equations. I underestimated just how complicated they would be to
    % implement, even with just that one extra constraint, so was trying to get
    % to grips with how they would link into each other in code. It turned out
    % actually to be useful/integral to write code which didn't include lots of
    % duplication and/or inefficiently passing way too many variables around.
    %
    % The `m2b' argument determines which branch of the algorithm the function
    % is being called in. If m2b = true, it's in move_to_bound, otherwise it's
    % in becomes_free. This changes how the b vector is calculated.
    %
    % See also, ???  % TODO lookup what MATLAB standard is when this is empty 

    % set default value for tolerance
    if (nargin < 13); tol = 1e-10; end  % TODO check this is sensible

    % want D and S to be indexes in free assets, not indexes in all the assets
    F_D = subindex(D, F);
    F_S = subindex(S, F);

    % matrix entries on the LHS of the multipler linear system
    a = sum(sum(invcovarF(F_D,F_S)));
    b = sum(sum(invcovarF(F_D,F_D)));
    c = sum(sum(invcovarF(F_S,F_S)));
    % d = a, so doesn't need repeating

    % determinant of the multiplier linear system
    determinant = a^2 - b*c;
    % BUG work out when this occurs; could be impossible in well formed problems
    if (abs(determinant) < tol)
        % TODO work out how to handle situations
        disp("ERROR! Linear system has no solutions")
    end

    % column sums of the four quadrants of invcovarF
    eF_DinvcovarF_DD = sum(invcovarF(F_D, F_D), 1);
    eF_DinvcovarF_DS = sum(invcovarF(F_D, F_S), 1);
    eF_SinvcovarF_SD = sum(invcovarF(F_S, F_D), 1);
    eF_SinvcovarF_SS = sum(invcovarF(F_S, F_S), 1);

    % TODO refactor to involve code duplication
    % components of x & y are reused elsewhere so split up
    x_p1 = 0.5;
    y_p1 = 0.5;
    % if B empty we get massive cancellations, so only consider non-empty case
    if ~isempty(B)
        % for x & y also need sub-indexed wB
        B_D = subindex(D, B);
        B_S = subindex(S, B);
        wB = w(B);
        wB_S = wB(B_S);
        wB_D = wB(B_D);
        % first part of x and y when B non-empty
        if isempty(B_D)
            printf("MU: B_D is empty, B_S is non-empty\n")
            % cancellations if no bounded dams
            x_p1 = x_p1 + eF_DinvcovarF_DD*(covarFB(F_D, B_S)*wB_S);
            y_p1 = y_p1 - sum(wB_S, 1) + eF_SinvcovarF_SS*(covarFB(F_S, B_S)*wB_S);
        elseif isempty(B_S)
            printf("MU: B_S is empty, B_D is non-empty\n")
            % cancellations if no bounded sires
            x_p1 = x_p1 - sum(wB_D, 1) + eF_DinvcovarF_DD*(covarFB(F_D, B_D)*wB_D);
            y_p1 = y_p1 + eF_SinvcovarF_SS*(covarFB(F_S, B_D)*wB_D);
        else
            printf("MU: B is non-empty\n")
            % no cancellations if we've bounded dams and sires
            x_p1 = x_p1 - sum(wB_D, 1) + eF_DinvcovarF_DD*(covarFB(F_D, B_S)*wB_S) + eF_DinvcovarF_DD*(covarFB(F_D, B_D)*wB_D);
            y_p1 = y_p1 - sum(wB_S, 1) + eF_SinvcovarF_SS*(covarFB(F_S, B_S)*wB_S) + eF_SinvcovarF_SS*(covarFB(F_S, B_D)*wB_D);
        end
    end
    printf("MU: B is empty\n")

    x_p2 = eF_DinvcovarF_DS*muF(F_S) + eF_DinvcovarF_DD*muF(F_D);
    y_p2 = eF_SinvcovarF_SS*muF(F_S) + eF_SinvcovarF_SD*muF(F_D);

    % final multiplier update is much simpler
    if (lam_current == 0)
        printf("MU: lambda is zero\n")
        % already know determ
        gam = (a*x_p1 - b*y_p1)/determinant;
        del = (a*y_p1 - c*x_p1)/determinant;
        lambda = 0;
        C = NaN;  % not needed
        return
    end
    printf("MU: lambda is non-zero\n")

    % derivatives of previous multipliers with respect to lambda
    dGam_dLam = (b*y_p2 - a*x_p2)/determinant;
    dDel_dLam = (c*x_p2 - a*y_p2)/determinant;

    % entry i of C gives whether free asset i is moving toward its upper or lower bound
    C = invcovarF * ([repmat(dGam_dLam, size(F_S)), repmat(dDel_dLam, size(F_D))] + muF);

    % lambda update equation is separated into chunks for readability. the first
    % of these chunks is determined by whether we're considering assets becoming free
    % or moving to a bound.

    if m2b
        bound = zeros(length(w), 1);  % b vector

        % in move_to_bound, so b(i) determined based on derivative
        for j = 1:length(F)
            % need to index outer matrix, as per NOTE A1
            i = F(j);
            if C(j) > 0
                % asset moving towards its upper bound
                bound(i) = ub(i);
            elseif C(j) < 0
                bound(i) = lb(i);
            else  % C(j) == 0
                % since more than one free variable, C(j) == 0 iff all mu(i) equal
                continue
            end
        end
        % precalculate lam_num_p1
        if isempty(B)
            printf("MU: MB: B is empty\n")
            lam_num_p1 = determinant*bound(F);
        else
            printf("MU: MB: B is non-empty\n")
            lam_num_p1 = determinant*(bound(F) + invcovarF*covarFB*w(B));
        end
    else
        % in becomes_free, so b(i) = w(i)
        if isempty(B)
            printf("MU: BF: B is empty\n")
            lam_num_p1 = determinant*w(F);
        else
            printf("MU: BF: B is non-empty\n")
            lam_num_p1 = determinant*(w(F) + invcovarF*covarFB*w(B));
        end
    end

    lam_num_p2 = [repmat(b*y_p1 - a*x_p1, size(F_S)), repmat(c*x_p1 - a*y_p1, size(F_D))]';
    lam_den_p1 = [repmat(b*y_p2 - a*x_p2, size(F_S)), repmat(c*x_p2 - a*y_p2, size(F_D))]';

    lam_den_p2 = invcovarF*(lam_den_p1 + determinant*muF);
    lambda = (lam_num_p1 + invcovarF*lam_num_p2)./lam_den_p2;

    % matrix entries on the RHS of the multipler linear system
    x = x_p1 - lambda*x_p2;
    y = y_p1 - lambda*y_p2;

    % update multipliers by direct calculateion of linear system solutions
    gam = (a*x - b*y)/determinant;
    del = (a*y - c*x)/determinant;

    % NOTE addityional x, y, gam, and del calculations are being done here
    % than necessary (because lambda is a vector and we only choose one).
    % Would be worth checking if we can make this better
end
