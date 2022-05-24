% +----------------------------------------------------------------------+
% | BELOW IS SCRAP AREA TESTING AREA, THESE FUNCTIONS ARE NOT FUNCTIONAL |
% | AND JUST HERE WHILE I'M WORKING ON MY CRAP, TRYING TO GET IT TO WORK |
% +----------------------------------------------------------------------+

function [gam, del, C, lambda] = multiplier_update(F, B, S, D, w, invcovarF, covarFB, muF, lam_current, tol=1e-10)
    % MULTIPLIER_UPDATE handle the updating of lagrangian multipliers
    % 
    % TODO work out if this is actually needed or can be simplified away
    % This was mostly me just working through the full set of multiplier update
    % equations trying to get to grips with how they would link into each other
    % in code. I'm not even sure it will actually be functional!
    %
    % See also, ???

    % want D and S to be indexes in free assets, not indexes in all the assets
    F_D = subindex(D, F);
    F_S = subindex(S, F);

    % matrix entries on the LHS of the multipler linear system
    a = sum(sum(invcovarF(F_D,F_S)));
    b = sum(sum(invcovarF(F_D,F_D)));
    c = sum(sum(invcovarF(F_S,F_S)));
    % d = a, so doesn't need repeating

    % determinant of the multiplier linear system
    determinant = a**2 - b*c;
    % TODO work out when this occurs; could be impossible in well formed problems
    if (abs(determinant) < tol)
        % TODO work out how to handle situations
        printf("ERROR! Linear system has no solutions")
    end

    % column sums of the four quadrants of invcovarF
    eF_DinvcovarF_DD = sum(invcovarF(F_D, F_D), 1);
    eF_DinvcovarF_DS = sum(invcovarF(F_D, F_S), 1);
    eF_SinvcovarF_SD = sum(invcovarF(F_S, F_D), 1);
    eF_SinvcovarF_SS = sum(invcovarF(F_S, F_S), 1);

    % components of x & y are reused elsewhere so split up
    if isempty(B)
        % if B empty we get massive cancellations
        x_p1 = 0.5;
        y_p1 = 0.5;
    else
        % for x & y also need sub-indexed wB
        wB_S = w(B)(subindex(S, B));
        wB_D = w(B)(subindex(D, B));
        % first part of x and y when B non-empty
        x_p1 = 0.5 - sum(wB_D, 1) + eF_DinvcovarF_DD*(covarFB(F_D, F_S)*wB_S) + eF_DinvcovarF_DD*(covarFB(F_D, F_D)*wB_D);
        y_p1 = 0.5 - sum(wB_S, 1) + eF_SinvcovarF_SS*(covarFB(F_S, F_S)*wB_S) + eF_SinvcovarF_SS*(covarFB(F_S, F_D)*wB_D);
    end

    x_p2 = eF_DinvcovarF_DS*muF(F_S) + eF_DinvcovarF_DD*muF(F_D);
    y_p2 = eF_SinvcovarF_SS*muF(F_S) + eF_SinvcovarF_SD*muF(F_D);

    % matrix entries on the RHS of the multipler linear system
    x = x_p1 - lam_current*x_p2;
    y = y_p1 - lam_current*y_p2;

    % calculate solutions of linear system directly; verified det != 0 already.
    gam = (a*x - b*y)/determinant;
    del = (a*y - c*x)/determinant;

    % NOTE calculating gam, del, x, or y is *not* a pre-requisite for C or
    % lambda. All we need for those are a, b, c, d and {x,y}_p{1,2}. It's not
    % a massive computational saving, but it's something.

    % derivatives of those updated multipliers with respect to lambda
    dGam_dLam = (b*y_p2 - a*x_p2)/determinant;
    dDel_dLam = (c*x_p2 - a*y_p2)/determinant;

    % entry i of C gives whether free asset i is moving toward its upper or lower bound 
    C = invcovarF * ([repmat(dGam_dLam, size(S)); repmat(dDel_dLam, size(D))] + muF);

    % lambda update equation is separated into chunk. the first of these chunks
    % is determined by whether we're considering assets becoming free or moving
    % to a bound.

    % BEGIN BLOCK
    % if in moves_to_bound, b(i) determined based on derivative
    b = zeros(size(w));
    for i = B
        if C(i) > 0
            % asset moving towards its upper bound
            b(i) = ub(i);
        elseif C(i) < 0
            b(i) = lb(i);
        else  % C(i) == 0
            % since more than one free variable, C(j) == 0 iff all mu(i) equal
            continue
        end
    end
    if isempty(B)
        lam_num_p1 = determinant*b;
    else
        lam_num_p1 = determinant*(b + invcovarF*covarFB*w(B));
    end
    % END BLOCK

    % BEGIN BLOCK
    % if in becomes_free function, b(i) = w(i)
    if isempty(B)
        lam_num_p1 = determinant*w;
    else
        lam_num_p1 = determinant*(w + invcovarF*covarFB*w(B));
    end
    % END BLOCK

    lam_num_p2 = [repmat(b*y_p1 - a*x_p1, size(S)); repmat(c*x_p1 - a*y_p1, size(D))];
    lam_den_p1 = [repmat(b*y_p2 - a*x_p2, size(S)); repmat(c*x_p2 - a*y_p2, size(D))];
    lam_den_p2 = invcovarF*(lam_den_p1 + determinant*muF);
    lam_new = (lam_num_p1 + lam_num_p2)/lam_den_p2;  % Needs to be vectorised
end


function [dGam_dLam, dDel_dLam, C] = derivative_update(F, B, S, D, w, invcovarF, covarFB, muF, tol=1e-10)
    % DERIVATIVE_UPDATE handle the CLA case where an asset becomes free
    %
    % Calculating the derivative is complicated in the two-half-constraints
    % version of the algorithm. Given it's done multiple times, this function
    % is included to handle the heavy lifting here.
    %
    % See also, MULTIPLIER_UPDATE

    % want D and S to be indexes in free assets, not indexes in all the assets
    F_D = subindex(D, F);
    F_S = subindex(S, F);

    % matrix entries on the LHS of the multipler linear system
    a = sum(sum(invcovarF(F_D,F_S)));
    b = sum(sum(invcovarF(F_D,F_D)));
    c = sum(sum(invcovarF(F_S,F_S)));
    % d = a, so doesn't need repeating

    % determinant of the multiplier linear system
    determinant = a**2 - b*c;
    % TODO work out when this occurs; could be impossible in well formed problems
    if (abs(determinant) < tol)
        % TODO work out how to handle situations
        printf("ERROR! Linear system has no solutions")
    end

    % column sums of the four quadrants of invcovarF
    eF_DinvcovarF_DD = sum(invcovarF(F_D, F_D), 1);
    eF_DinvcovarF_DS = sum(invcovarF(F_D, F_S), 1);
    eF_SinvcovarF_SD = sum(invcovarF(F_S, F_D), 1);
    eF_SinvcovarF_SS = sum(invcovarF(F_S, F_S), 1);

    % components of x & y are reused elsewhere so split up
    if isempty(B)
        % if B empty we get massive cancellations
        x_p1 = 0.5;
        y_p1 = 0.5;
    else
        % for x & y also need sub-indexed wB
        wB_S = w(B)(subindex(S, B));
        wB_D = w(B)(subindex(D, B));
        % first part of x and y when B non-empty
        x_p1 = 0.5 - sum(wB_D, 1) + eF_DinvcovarF_DD*(covarFB(F_D, F_S)*wB_S) + eF_DinvcovarF_DD*(covarFB(F_D, F_D)*wB_D);
        y_p1 = 0.5 - sum(wB_S, 1) + eF_SinvcovarF_SS*(covarFB(F_S, F_S)*wB_S) + eF_SinvcovarF_SS*(covarFB(F_S, F_D)*wB_D);
    end

    x_p2 = eF_DinvcovarF_DS*muF(F_S) + eF_DinvcovarF_DD*muF(F_D);
    y_p2 = eF_SinvcovarF_SS*muF(F_S) + eF_SinvcovarF_SD*muF(F_D);

    % Calculating gam, del, x, or y is *not* a pre-requisite for C so they're
    % ommitted here to for efficincy. It's not a huge saving but it's something.

    % derivatives of the new multipliers with respect to lambda
    dGam_dLam = (b*y_p2 - a*x_p2)/determinant;
    dDel_dLam = (c*x_p2 - a*y_p2)/determinant;

    % entry i of C gives whether free asset i is moving toward its upper or lower bound 
    C = invcovarF * ([repmat(dGam_dLam, size(S)); repmat(dDel_dLam, size(D))] + muF);
end
