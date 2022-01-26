version = '0.3.1';
printf(['OctCLA v' version '\n'])

% This file octcla.m is part of OctCLA, Copyright (c) 2022 Josh Fogg, released
% under the MIT License. For more info see https://github.com/Foggalong/OctCLA


% UTILITY FUNCTIONS
% =================
% These functions aren't actually part of the algorithm themselves, but they
% carry out calculations or operations the results of which the algorithm uses.

function [imax, xmax] = argmax(x, a=inf, tol=1e-10)  % TODO check this is sensible
    % ARGMAX return index and value of highest item in a list
    %
    % Takes a vector x as input and then returns the index imax and value xmax
    % of the highest value in x, or NA if x is empty. Can also take an optional
    % second input a, restricting the chosen x(j) to satisfy x(j) < a.
    %
    % See also, MAX.

    % starting values catch the case where x empty
    imax = 0;
    xmax = -inf;
    % TODO see if this function can be vectorised
    for i = 1:length(x)
        % if a=inf second check always satified
        if (tol < x(i)-xmax && tol < a-x(i))
            xmax = x(i); imax = i;
        end
    end
    % if no max found, return NA
    if imax == 0
        imax = NA; xmax = NA;
    end
end

function newInv = inverse_grow(invA, a, alpha)
    % INVERSE_SHRINK adjust the inverse if gaining a row and column
    %
    % Takes the inverse of a symmetric matrix, a column vector, and a scalar
    % as inputs and then returns the inverse of the original matrix if its
    % dimensions were both increased by 1, i.e. if newA = [A, a; a, alpha]
    % it takes A^-1, a, and alpha as input and returns newA^-1.
    %
    % See also, CALCULATE_TURNINGPOINTS

    % shortcut variabels
    c = invA*a;
    beta = 1/(alpha - c'*a);
    % preallocate the new inverse
    newInv = zeros(size(invA)+1);
    % update the top left block
    newInv(1:end-1, 1:end-1) = invA + beta*(c*c');
    % update the top right block
    newInv(1:end-1, end) = -beta*c;
    % update the bottom left block
    newInv(end, 1:end-1) = -beta*c';
    % update the bottom right block
    newInv(end, end) = beta;
end

function newInv = inverse_shrink(invA, i)
    % INVERSE_SHRINK adjust the inverse if removing row and column i
    %
    % Takes the inverse of a symmetric matrix and an index i as input and then
    % returns the inverse of the original matrix if its dimensions were both
    % decreased by 1 through removal of row and column i. That is, if after
    % permuting the ith columns and rows to the end, A = [newA, a; a, alpha]
    % then the function takes A^-1 (and i) as input and returns newA^-1.
    %
    % See also, INVERSE_GROW

    % permute rows and cols of invA so that ith row and col are last
    invAperm = invA([1:(i-1) (i+1):end i], [1:(i-1) (i+1):end i]);
    % shortcut variables
    B = invAperm(1:end-1, 1:end-1);
    b = invAperm(1:end-1, end);
    beta = invAperm(end, end);
    % calculate new inverse
    newInv = B - (b*b')/beta;
end


% KKT CHECK FUNCTIONS
% ===================
% These two functions carry out a KKT check on a solution to a problem with
% the corresponding input variables. They make use of a vector d calculated
% by the ALGORITHM FUNCTIONS which contains the value of the bound a given
% weight is moving towards (itself calculated using the derivative).

function result = kkt_full(lb, ub, d, w_t, covar, lam, mu, gam,
                           debug=false, tol=1e-10)  % TODO check this is sensible
    % KKT_FULL perform a full KKT check on a constrained CLA problem
    %
    % Takes lambda (lam), gamma (gam), and derivative vector (d) values
    % associated with a solution (w_t) to a portfolio optimisation problem
    % with given covariance matrix (covar), expected return vector (mu), and
    % bound constraints (lb, ub) and then returns true if that solution passes
    % the KKT conditions and false otherwise.
    %
    % Takes an optional bool variable `debug` which prints which condition was
    % broken and the relevant values. Takes another optional float input which
    % is the tollerance to which the conditions are verified.
    %
    % See also, KKT_EQUIV

    % Check KKT(2.2): weights sum to one
    if sum(w_t) - 1 > tol
        if debug
            w_t = w_t
            printf('Broke KKT(2.2)\n');
        end
        result = false; return
    end

    % The ith entry of vector d contains the value of the bound which the ith
    % weight was moving toward in the solution. Thus the diagonal matrix C is
    % such that the (i,i)th entry is +1 if the ith weight was moving toward the
    % lower bound and -1 if it was moving toward the upper bound.
    C = diag([d==lb] - [d==ub]);

    % Check KKT(3): weights within bounds
    if C*w_t-d > tol
        if debug
            d_Cw = d-C*w_t
            printf('Broke KKT(3)\n');
        end
        result = false; return
    end

    % KKT(1) is assumed through this statement
    tau = C*(covar*w_t - lam*mu - gam*ones(size(mu)));

    % Check KKT(4): Tau positive
    if tau < -tol
        if debug
            tau = tau
            printf('Broke KKT(4)\n');
        end
        result = false; return
    end

    % Check KKT(5): tau'*(Cw-d) = 0
    if abs(tau'*(C*w_t - d)) > tol
        if debug
            d = d
            Cw_d = C*w_t-d
            tau = tau
            tCw_d = tau'*(C*w_t - d)
            printf('Broke KKT(5)\n');
        end
        result = false; return
    end

    result = true;
end

function result = kkt_equiv(covar, w, lam, gam, mu, B, F,
                           debug=false, tol=1e-10)  % TODO check this is sensible
    % KKT_EQUIV unconstrained-equivalent KKT check on a constrained CLA problem
    %
    % Takes a lambda (lam), a gamma (gam), and (F)ree and (B)ounded weight lists
    % associated with a solution (w) to an unconstrained portfolio optimisation
    % problem, which is the equivalent problem to the constrained problem of
    % KKT_FULL, and then returns true if that solution passes the KKT conditions
    % and false otherwise.
    %
    % Takes an optional bool variable `debug` which prints which condition was
    % broken and the relevant values. Takes another optional float input which
    % is the tollerance to which the conditions are verified.
    %
    % See also, KKT_FULL

    % Check KKT(2.2): weights sum to one
    if sum(w) - 1 > tol
        if debug
            w
            printf('Broke KKT(2.2)\n');
        end
        result = false; return
    end

    % Check KKT(1): Lagrangian is zero
    onesF = ones(size(mu(F)));
    delwL = covar(F,F)*w(F) + covar(F,B)*w(B) - lam*mu(F) - gam*onesF;
    if abs(delwL) > tol
        if debug
            delwL
            printf('Broke KKT(1)\n');
        end
        result = false; return
    end

    result = true;
end


% ALGORITHM FUNCTIONS
% ===================
% These four functions make up the key part of the algorithm; finding the
% starting solution, handling the case when an asset moves to its bound,
% handling the case when an asset become free, and then calculating the
% turning points themselves through CLA.

function [F, B, w] = starting_solution(mu, lb, ub)
    % STARTING_SOLUTION return starting solution for CLA
    %
    % Takes a vector mu of expected returns, a vector lb of lower bounds on
    % assest weights, and a vector ub of upper bounds on asset weights, then
    % returns the starting solution for CLA in the form of an index list F of
    % free assests and starting weight vector w.
    %
    % See also, CALCULATE_TURNINGPOINTS

    % start with all assets on their lower bound
    w = lb;
    % increase assest weights in descending order of expected return
    i = argmax(mu);
    while sum(w) < 1
        i_free = i;
        w(i) = min(ub(i), lb(i)+1-sum(w));
        i = argmax(mu, mu(i));
    end
    % only one asset starts free
    F = [i_free];
    % all other assets start on bounds
    B = 1:length(mu);
    B = B(B ~= i_free);
end

% NOTE added a B argument which isn't in the original, need to check if it can be removed
% NOTE also added a d vector which is related to the KKD conditions check
function [ins, lam_ins, b_ins, d] = move_to_bound(mu, covar, invcovarF, lb, ub, F, B, lam_current, w)
    % MOVE_TO_BOUND handle the CLA case where an asset moves to its bound
    %
    % Takes a vector mu of expected returns, a covariance matrix covar, a
    % vector lb of lower bounds on assest weights, and a vector ub of upper
    % bounds on asset weights, then returns the i, lambda, and bound of the
    % asset which would move to its bound.
    %
    % See also, CALCULATE_TURNINGPOINTS

    % a sole free asset cannot move to bound
    if length(F) == 1
        ins = b_ins = d = NA; lam_ins = -inf;
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
            % since more than one free variable, C(j) == 0 iff all mu(i) equal
            continue  % TODO check how properly to handle this
        end

        % calculate lambda using shortcuts TODO this could really be tidied up
        if isempty(B)
            lami_p2 = sum(sum(invcovarF))*b(i);
        else
            lami_p2 = sum(sum(invcovarF))*(b(i)+(invcovarF*(covarFB*w(B)))(j));
        end

        lam(i) = (lam_p1(j)-lami_p2)/C(j);
    end
    d = b;
    d(B) = w(B);

    [ins, lam_ins] = argmax(lam, lam_current);
    if isnan(ins)
        b_ins = NA;
    else
        b_ins = b(ins);
    end
end

% NOTE added a B argument which isn't in the original, need to check if it can be removed
function [outs, lam_outs, d] = becomes_free(mu, covar, invcovarF, lb, ub, F, B, lam_current, w)
    % BECOMES_FREE handle the CLA case where an asset becomes free
    %
    % Takes a vector mu of expected returns, a covariance matrix covar, a
    % vector lb of lower bounds on assest weights, and a vector ub of upper
    % bounds on asset weights, then returns the i and lambda of the asset
    % which moves away from its bound (i.e. which becomes free).
    %
    % See also, CALCULATE_TURNINGPOINTS

    % skip proceedure if all assets are free
    if (length(F) == length(mu))
        outs = d = NA; lam_outs = -inf;
        return
    end

    D   = zeros(length(mu), length(mu));  % matrix of potential d vectors
    lam = zeros(length(mu), 1);           % lambda vector

    for i = B
        % update the inverse
        a = covar(F, i);
        alpha = covar(i, i);
        invcovarFi = inverse_grow(invcovarF, a, alpha);
        % free weight i
        Fi = [F, i];     % F = Fu{i}
        Bi = B(B ~= i);  % B = B\{i}
        % need to index outer matrix, as per NOTE A1
        j = length(Fi);  % Fi[j] = i; i last element in Fi by construction

        % calculate derivative TODO find a way to make this somewhat nicer
        Ci1 = -sum(sum(invcovarFi))*(invcovarFi*mu(Fi));
        Ci2 = (sum(invcovarFi)*mu(Fi))*sum(invcovarFi, 2);
        Ci = Ci1 + Ci2;

        % saving d vector for KKT conditions
        for l = 1:length(Fi)
            k = Fi(l);
            if Ci(l) > 0; D(k, i) = ub(l); end
            if Ci(l) < 0; D(k, i) = lb(l); end
        end
        D(Bi, i) = w(Bi);

        % handle case in NOTE A2
        if isempty(Bi)
            lami_p1 = sum(invcovarFi, 2)(j);
            lami_p2 = sum(sum(invcovarFi))*w(i);
        else
            % calculate lambda using shortcuts TODO this could really be tidied up
            lami_p1 = (1-sum(w(Bi))+sum(invcovarFi)*(covar(Fi,Bi)*w(Bi)))*(sum(invcovarFi,2)(j));
            lami_p2_q2 = invcovarFi*(covar(Fi,Bi)*w(Bi));
            lami_p2 = sum(sum(invcovarFi))*(w(i)+lami_p2_q2(j));
        end

        lam(i) = (lami_p1-lami_p2)/Ci(j);
    end

    [outs, lam_outs] = argmax(lam, lam_current);
    if outs ~= NA
        d = D(:, outs);
    else
        d = NA;
    end
end

function ws = calculate_turningpoints(mu, covar, lb, ub, KKT=1)
    % CALCULATE_TURNINGPOINTS return portfolio weights at CLA turning points
    %
    % Takes a vector mu of expected returns, a covariance matrix covar, a
    % vector lb of lower bounds on assest weights, and a vector ub of upper
    % bounds on asset weights, then calculates the turning points of the
    % corresponding portfolio using the Critical Line Algorithm. The function
    % returns a matrix of the portfolio asset weights at each of those points.
    %
    % Takes an optional integer argument `KKT` which indicates whether calculated
    % turning points should be verified through a KKT conditions check and if so
    % by what type of check (0: none, 1: simple, 2: full, 3: both). Default is 1.
    %
    % See also, STARTING_SOLUTION, BECOMES_FREE, MOVE_TO_BOUND

    % calculate starting solution
    [F, B, ws] = starting_solution(mu, lb, ub);

    % initial inversion, the only one calculated without shortcuts
    invcovarF = inv(covar(F,F));

    lam_current = inf;
    t = 1;  % so current w will be ws(:,t)

    while true
        % case a where a free asset moves to its bound
        [i_ins, lam_ins, b, d_ins] = move_to_bound(mu, covar, invcovarF, lb, ub, F, B, lam_current, ws(:,t));
        % case b where an asset on its bound becomes free
        [i_outs, lam_outs, d_outs] = becomes_free(mu, covar, invcovarF, lb, ub, F, B, lam_current, ws(:,t));

        if (i_ins ~= NA || i_outs ~= NA)
            lam_current = max(lam_ins, lam_outs);
            % if lam < 0 the risk is increasing again, make lam = 0 the last iteration
            if lam_current < 0; lam_current = 0; end

            % add an additional column to ws for the new asset weights
            ws = [ws, zeros(length(mu),1)];
            ws(B, t+1) = ws(B, t);
            t = t+1;

            % update gamma TODO make this neater
            if isempty(B)
                gam_top = -lam_current*sum(invcovarF)*mu(F) + 1;
            else
                gam_top = -lam_current*sum(invcovarF)*mu(F) + 1 - sum(ws(B,t)) + sum(invcovarF)*(covar(F,B)*ws(B,t));
            end
            gam = gam_top/sum(sum(invcovarF));

            % update w_F^(t) TODO make this neater
            if isempty(B)
                ws(F,t) = gam*sum(invcovarF,2) + lam_current*(invcovarF*mu(F));
            else
                ws(F,t) = -invcovarF*(covar(F,B)*ws(B,t)) + gam*sum(invcovarF,2) + lam_current*(invcovarF*mu(F));
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
                ws(i_ins, t) = b;   % w_i_inside^(t) = b TODO is this line needed?
                d = d_ins;
            else
                % update the inverse
                a = covar(F, i_outs);
                alpha = covar(i_outs, i_outs);
                invcovarF = inverse_grow(invcovarF, a, alpha);
                % free weight i_outs
                F = [F, i_outs];     % F = Fu{i}
                B = B(B ~= i_outs);  % B = B\{i}
                d = d_outs;
            end

            % KKT CHECKS
            % NOTE doesn't check lam=0 soln since it's not (necessarily) a TP.
            if (KKT > 0)
                printf('Checking KKT conditions...')
                check = 0;

                if (KKT == 1) || (KKT == 3) 
                    check += kkt_full(lb, ub, d, ws(:, t), covar, lam_current, mu, gam, debug=true);
                end

                if (KKT == 2) || (KKT == 3)
                    check += kkt_equiv(covar, ws(:,t), lam_current, gam, mu, B, F, debug=true);
                end

                if (KKT == 3 && check == 2) || (KKT ~= 3 && check == 1)
                    printf(' passed!\n')
                else
                    if (KKT ~= 3) || (check ~= 1)
                        printf(' checks failed!\')
                    else
                        printf(' uh oh, only one check failed!\n')
                    end
                    % if DEBUG
                    w = ws(:,t)'
                    all_w = ws'
                    cond(covar)
                    exit(1)
                    % end DEBUG
                end
            end

        else
            % if i_ins and i_outs are NA then we're done
            break
        end
    end
    % only return w2, w3, etc since w0 and w1 coincide TODO work out why
    ws = ws(:, 2:end)';
end
