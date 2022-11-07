disp('OctCLA Genetics v0.4.1')

% This file octcla-gen.m is part of OctCLA, Copyright (c) 2022 Josh Fogg, released
% under the MIT License. For more info see https://github.com/Foggalong/OctCLA


% UTILITY FUNCTIONS
% ==================
% These functions aren't actually part of the algorithm themselves, but they
% carry out calculations or operations the results of which the algorithm uses.

% This file depends on the following functions from the standard implementation:
%  - MAX_BOUNDED maximum value in x satisfying x(i) < x_bound
%  - INVERSE_SHRINK adjust the inverse if gaining a row and column
%  - INVERSE_SHRINK adjust the inverse if removing row and column i
% source('../octcla.m')  % NOTE command only GNU Octave compatible

function [x_max, i_max] = max_conditional(x, x_bound, index_set, tol)
    % MAX_CONDITIONAL return value and index of highest item meeting conditions
    %
    % Takes a vector x, a number x_bound, and a set of indices index_set as inputs
    % and returns the value x_max and index i_max of the highest value in x
    % satisfying conditions that i_max is in index_set and x(i_max) < x_bound.
    % Note if index_set = 1:length(x) the function is equivalent to max_bounded.
    % Takes an optional value tol which specifies tolerance (default: 10^-10).
    %
    % See also, MAX, MAX_BOUNDED

    % set default value for tolerance
    if (nargin < 4); tol = 1e-10; end

    % starting values catch the case where x empty
    x_max = -inf;
    i_max = 0;

    % only interested in indicies from index_set so iterate over them
    for i = index_set
        % if x_bound=inf second check always satified
        if (tol < x(i)-x_max && tol < x_bound-x(i))
            x_max = x(i); i_max = i;
        end
    end

    % if no max found, return NaN
    if i_max == 0
        x_max = NaN; i_max = NaN;
    end
end

% HACK Sometimes we will use S and D to index the full set of assets and sometimes
% we will only be indexing the free or bounded assets. A hack to let us do that
% now is this function which lets us convert from one to the other.

function output = subindex(X, Y)
    % SUBINDEX returns locations of elements of X appearing in Y
    %
    % Given index vectors X and Y of some other vector V, return an index
    % vector of Y giving locations of elements of X that appear in Y.
    %
    % See also, ISMEMBER.

    [bool_XinY, ind_XinY] = ismember(X, Y);
    output = ind_XinY(bool_XinY);
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

function [F, B, w] = starting_solution_gen(mu, lb, ub, S, D, tol)
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

    % set default value for tolerance
    if (nargin < 6); tol = 1e-10; end  % TODO check this is sensible

    % start with all assets on their lower bound
    w = lb;

    % setting condition to inf ensures 1st value satisfies mi(i) < mu_max
    mu_max = inf;
    % increase sire weights in descending order of expected return
    while (abs(sum(w(S)) - 0.5) > tol)
        [mu_max, i] = max_conditional(mu, mu_max, S);
        w(i) = min(ub(i), lb(i)+0.5-sum(w(S)));
    end
    % only one sire starts free
    F = [i];
    B = S(S ~= i);

    % setting condition to inf ensures 1st value satisfies mi(i) < mu_max
    mu_max = inf;
    % increase dam weights in descending order of expected return
    while (abs(sum(w(D)) - 0.5) > tol)
        [mu_max, i] = max_conditional(mu, mu_max, D);
        w(i) = min(ub(i), lb(i)+0.5-sum(w(D)));
    end
    % only one dam starts free
    F = [F, i];
    % all other dams start on bounds
    B = [B, D(D ~= i)];
end


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
        % already know determ
        gam = (a*x_p1 - b*y_p1)/determinant;
        del = (a*y_p1 - c*x_p1)/determinant;
        lambda = 0;
        C = NaN;  % not needed
        return
    end

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
            lam_num_p1 = determinant*bound(F);
        else
            lam_num_p1 = determinant*(bound(F) + invcovarF*covarFB*w(B));
        end
    else
        % in becomes_free, so b(i) = w(i)
        if isempty(B)
            lam_num_p1 = determinant*w(F);
        else
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
    b   = zeros(length(mu), 1);  % holds which bounds being moved towards
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
    elseif ((length(F_D) == 1) && (ismember(ins, D)))
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

function ws = calculate_turningpoints_gen(mu, covar, lb, ub, S, D, KKT)
    % CALCULATE_TURNINGPOINTS_GEN return portfolio weights at genetics CLA turning points
    %
    % Takes a vector of expected returns (mu), a covariance matrix (covar), a
    % vector of lower bounds on assest weights (lb), a vector of upper bounds
    % on asset weights (ub), a set of sire indexes (S), and a set of dam indexes
    % (D) and then calculates the turning points of the corresponding genetics
    % portfolio using the Critical Line Algorithm. The function returns a matrix
    % of the portfolio asset weights at each of those points (ws).
    %
    % Takes an optional integer argument (KKT) which indicates whether calculated
    % turning points should be verified through a KKT conditions check and if so
    % by what type of check (0: none, 1: simple (default), 2: full, 3: both).
    %
    %
    % See also, STARTING_SOLUTION_GEN, BECOMES_FREE_GEN, MOVE_TO_BOUND_GEN

    % set default value for KKT
    if (nargin < 7); KKT = 1; end

    % calculate starting solution
    [F, B, ws] = starting_solution_gen(mu, lb, ub, S, D);

    % initial inversion, the only one calculated without shortcuts
    % TODO this is almost guaranteed to be a 2x2 matrix, should be able to do this better 
    invcovarF = inv(covar(F,F));

    lam_current = inf;
    t = 1;  % so current w will be ws(:,t)

    while true
        % case a where a free asset moves to its bound
        % DEBUG - uncomment once done 
        [i_ins, lam_ins, gam_ins, del_ins, b, d_ins] = move_to_bound_gen(mu, covar, invcovarF, lb, ub, F, B, S, D, lam_current, ws(:,t), KKT);
        % case b where an asset on its bound becomes free
        [i_outs, lam_outs, gam_outs, del_outs, d_outs] = becomes_free_gen(mu, covar, invcovarF, mu(F), lb, ub, F, B, S, D, lam_current, ws(:,t), KKT);

        if (i_ins ~= NaN || i_outs ~= NaN)
            lam_current = max(lam_ins, lam_outs);

            % if lam < 0 the risk is increasing again, make lam = 0 the last iteration
            if lam_current < 0;
                lam_current = 0;
                % since value of lambda change need to recalculate gamma and delta
                [gam, del, C_null, lam_null] = multiplier_update(F, B, S, D, ws(F,t), invcovarF, covar(F,B), mu(F), lam_current, lb, ub, false);
                % NOTE don't actually need C or lambda; which makes choice of m2b redudent
            else
                % need to know which lambda won
                if lam_current == lam_ins  % can do without tol comparison since come from max
                    disp("Going inside")
                    gam = gam_ins;
                    del = del_ins;
                else
                    disp("Going outside")
                    gam = gam_outs;
                    del = del_outs;
                end
            end

            % add an additional column to ws for the new asset weights
            ws = [ws, zeros(length(mu),1)];
            ws(B, t+1) = ws(B, t);
            t = t+1;

            % update w_F^(t)
            F_D = subindex(D, F);
            F_S = subindex(S, F);
            ws(F,t) = lam_current*(invcovarF*mu(F)) + [gam*sum(invcovarF(F_S,F_S), 2)', del*sum(invcovarF(F_D,F_D), 2)']';
            % have an extra term unless b is empty
            if ~isempty(B)
                ws(F,t) = ws(F,t) - invcovarF*(covar(F,B)*ws(B,t));
            end

            % if lambda = 0 then risk is increasing again, can terminate
            if lam_current <= 0; break; end

            % update free and bounded asset index sets
            if (lam_ins > lam_outs)
                % update the inverse
                % TODO could this find be replaced?
                j = find(F==i_ins);  % need index in inverse, not full matrix
                invcovarF = inverse_shrink(invcovarF, j);
                % bound weight i_ins
                F = F(F ~= i_ins);  % F = F\{i}
                B = [B, i_ins];     % B = Bu{i}
                % ws(i_ins, t) = b;   % w_i_inside^(t) = b  % TODO is this line needed?
                % only need to update d if doing full KKT check
                if (KKT == 1) || (KKT == 3); d = d_ins; end
            else
                % update the inverse
                a = covar(F, i_outs);
                alpha = covar(i_outs, i_outs);
                invcovarF = inverse_grow(invcovarF, a, alpha);
                % free weight i_outs
                F = [F, i_outs];     % F = Fu{i}
                B = B(B ~= i_outs);  % B = B\{i}
                % only need to update d if doing full KKT check
                if (KKT == 1) || (KKT == 3); d = d_outs; end
            end

            % KKT CHECKS
            % NOTE doesn't check lam=0 soln since it's not (necessarily) a TP.

            if (KKT > 0)
                disp('Checking KKT conditions...')
                % ws'
                errors = 0;

                if (KKT == 1) || (KKT == 3)
                    % TODO update once associated function is implemented
                    % errors += ~kkt_full_gen(...);
                end

                if (KKT == 2) || (KKT == 3)
                    % TODO update once associated function is implemented
                    % errors += ~kkt_equiv(...);
                end

                if (errors == 0)
                    disp('...passed!')
                else
                    if (KKT == 3) && (errors == 1) 
                        disp('...uh oh, only one set of check failed!')
                    else
                        disp('...checks failed!')
                    end

                    w = ws(:,t)'
                    all_w = ws'
                    cond(covar)
                    exit(1)
                end
            end
        else
            % if i_ins and i_outs are NaN then we're done
            break
        end
    end
    % only return w2, w3, etc since w0 and w1 coincide  % TODO work out why
    ws = ws(:, 2:end)';
end
