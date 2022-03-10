version = '0.4.1';
printf(['OctCLA Genetics v' version '\n'])

% This file octcla-gen.m is part of OctCLA, Copyright (c) 2022 Josh Fogg, released
% under the MIT License. For more info see https://github.com/Foggalong/OctCLA


% UTILITY FUNCTIONS
% ==================
% These functions aren't actually part of the algorithm themselves, but they
% carry out calculations or operations the results of which the algorithm uses.

% This file depends on the following functions from the standard implementation:
%  - INVERSE_SHRINK adjust the inverse if gaining a row and column
%  - INVERSE_SHRINK adjust the inverse if removing row and column i
source('octcla.m')  % NOTE must come before redefining 

function [imax, xmax] = argmax(x, a=inf, X=1:length(x), tol=1e-10)  % TODO check this is sensible
    % ARGMAX return index and value of highest item in a list
    %
    % Takes a vector x as input and then returns the index imax and value xmax
    % of the highest value in x, or NA if x is empty. Can also take an optional
    % second input a, restricting the chosen x(j) to satisfy x(j) < a, an
    % optional third argument X, an index set from which j must come, and an
    % optional fourth argument tol for controlling the comparison tolerance.
    %
    % See also, MAX.

    % starting values catch the case where x empty
    imax = 0;
    xmax = -inf;
    % TODO see if this function can be vectorised
    for i = X
        % if a=inf second check always satified
        if (tol < x(i)-xmax && tol < a-x(i) && ismember(i, X))
            xmax = x(i); imax = i;
        end
    end
    % if no max found, return NA
    if imax == 0
        imax = NA; xmax = NA;
    end
end





% KKT CHECK FUNCTIONS
% ===================
% These two functions carry out a KKT check on a solution to a problem with
% the corresponding input variables. They make use of a vector d calculated
% by the ALGORITHM FUNCTIONS which contains the value of the bound a given
% weight is moving towards (itself calculated using the derivative).

% TODO implement replacement function for KKT_FULL

% TODO implement replacement function for KKT_EQUIV


% ALGORITHM FUNCTIONS
% ===================
% These four functions make up the key part of the algorithm; finding the
% starting solution, handling the case when an asset moves to its bound,
% handling the case when an asset become free, and then calculating the
% turning points themselves through CLA.

function [F, B, w] = starting_solution_gen(mu, lb, ub, S, D, tol=1e-10)
    % STARTING_SOLUTION_GEN return starting solution for CLA genetics problems
    %
    % Takes a vector of expected returns (mu), a vector of lower bounds on assest
    % weights (lb), a vector of upper bounds on asset weights (ub), an index set
    % for sires (S), and an index set for dams (D) as input, then returns the
    % starting solution for CLA in the form of an index list of (F)ree and
    % (B)ounded assests and starting weight vector (w). Also takes an optional
    % argument for controlling the tolerance for comparisons (tol).
    %
    % See also, STARTING_SOLUTION, CALCULATE_TURNINGPOINTS_GEN

    % start with all assets on their lower bound
    w = lb;

    % increase sire weights in descending order of expected return
    i = argmax(mu, a=inf, X=S);
    while (abs(sum(w(S)) - 0.5) > tol)
        i_free = i;
        w(i) = min(ub(i), lb(i)+0.5-sum(w(S)));
        i = argmax(mu, a=mu(i), X=S);
    end
    % only one sire starts free
    F = [i_free];
    B = S(S ~= i_free);

    % increase dam weights in descending order of expected return
    i = argmax(mu, a=inf, X=D);
    while (abs(sum(w(D)) - 0.5) > tol)
        i_free = i;
        w(i) = min(ub(i), lb(i)+0.5-sum(w(D)));
        i = argmax(mu, a=mu(i), X=D);
    end
    % only one dam starts free
    F = [F, i_free];
    % all other dams start on bounds
    B = [B, D(D ~= i_free)];
end


% TODO this won't work because the indexes sets D and S are for invcovar, not invcovarF
% will need to make a function which converts D to indexes in Fi


function output = subindex(X, Y)
    % given index vectors X and Y of some other vector V, return an index
    % vector of Y giving locations of elements of X that appear in Y.
    [bool_XinY, ind_XinY] = ismember(X, Y);
    output = ind_XinY(bool_XinY);
end



function [gam, del, C] = multiplier_update(F, B, S, D, invcovarF, covarFB, muF, lam_current, tol=1e-10)
    tol = 1e-10;
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

    % for x & y also need sub-indexed wB
    wB_S = w(B)(subindex(S, B));
    wB_D = w(B)(subindex(D, B));

    % column sums of the four quadrants of invcovarF
    eF_DinvcovarF_DD = sum(invcovarF(F_D, F_D), 1);
    eF_DinvcovarF_DS = sum(invcovarF(F_D, F_S), 1);
    eF_SinvcovarF_SD = sum(invcovarF(F_S, F_D), 1);
    eF_SinvcovarF_SS = sum(invcovarF(F_S, F_S), 1);

    % components of x & y are reused elsewhere so split up
    x_p1 = 0.5 - sum(wB_D, 1) + eF_DinvcovarF_DD*(covarFB(F_D, F_S)*wB_S) + eF_DinvcovarF_DD*(covarFB(F_D, F_D)*wB_D);
    y_p1 = 0.5 - sum(wB_S, 1) + eF_SinvcovarF_SS*(covarFB(F_S, F_S)*wB_S) + eF_SinvcovarF_SS*(covarFB(F_S, F_D)*wB_D);
    x_p2 = eF_DinvcovarF_DS*muF(F_S) + eF_DinvcovarF_DD*muF(F_D);
    y_p2 = eF_SinvcovarF_SS*muF(F_S) + eF_SinvcovarF_SD*muF(F_D);

    % matrix entries on the RHS of the multipler linear system
    x = x_p1 - lam_current*x_p2;
    y = y_p1 - lam_current*y_p2;

    % calculate solutions of linear system directly; verified det != 0 already.
    gam = (a*x - b*y)/determinant;
    del = (a*y - c*x)/determinant;

    % derivatives of those updated multipliers with respect to lambda
    dGam_dLam = (b*y_p2 - a*x_p2)/determinant;
    dDel_dLam = (c*x_p2 - a*y_p2)/determinant;

    % entry i of C gives whether free asset i is moving toward its upper or lower bound 
    C = invcovarF * ([repmat(dGam_dLam, size(S)); repmat(dDel_dLam, size(D))] + muF);

    lam_current = determinant* + invcovarF*covarFB*w(B))
end


mu = [0.2; 0.8; 0.4; 0.3];
lb = [0.2; 0.2; 0.2; 0.2];
ub = [0.5; 0.5; 0.5; 0.5];
S  = [1, 3];
D  = [2, 4];

starting_solution(mu, lb, ub, S, D);
