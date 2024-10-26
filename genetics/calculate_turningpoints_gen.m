% This file is part of OctCLA, Copyright (c) 2022 Josh Fogg, released
% under the MIT License. See: https://github.com/Foggalong/OctCLA

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
        [i_ins, lam_ins, gam_ins, del_ins, bound, d_ins] = move_to_bound_gen(mu, covar, invcovarF, lb, ub, F, B, S, D, lam_current, ws(:,t), KKT);
        % case b where an asset on its bound becomes free
        [i_outs, lam_outs, gam_outs, del_outs, d_outs] = becomes_free_gen(mu, covar, invcovarF, mu(F), lb, ub, F, B, S, D, lam_current, ws(:,t), KKT);

        if (~isnan(i_ins) || ~isnan(i_outs))
            lam_current = max(lam_ins, lam_outs);

            % if lam < 0 the risk is increasing again
            if lam_current < 0;
                % make lam = 0 the last iteration
                lam_current = 0;
                % since lambda changed, recalculate gamma and delta
                [gam, del, C_null, lam_null, b_null] = multiplier_update(F, B, S, D, ws(F,t), invcovarF, covar(F,B), mu(F), lam_current, lb, ub, false);
                % don't need C or lambda, so choice of m2b redudent
            else
                % can compare without tol since came from max
                if lam_current == lam_ins
                    gam = gam_ins; del = del_ins;    % disp("Going inside")
                else
                    gam = gam_outs; del = del_outs;  % disp("Going outside")
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
            % have an extra term unless B is empty
            if ~isempty(B)
                ws(F,t) = ws(F,t) - invcovarF*(covar(F,B)*ws(B,t));
            end

            % if lambda = 0, risk is increasing again: terminate
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
                % TODO is the below line needed?
                ws(i_ins, t) = bound;
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
                errors = 0;

                if (KKT == 1) || (KKT == 3)
                    errors = errors + ~kkt_full_gen(lb, ub, d, ws(:,t), covar, lam_current, mu, gam, del, S, D);
                end

                if (KKT == 2) || (KKT == 3)
                    errors = errors + ~kkt_partitioned_gen(covar, ws(:,t), lam_current, gam, del, mu, B, F, S, D);
                end

                if (errors == 0)
                    disp('...passed!')
                else
                    if (KKT == 3) && (errors == 1) 
                        disp('...uh oh, only one set failed!')
                    else
                        disp('...checks failed!')
                    end

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
    % only return w2, w3, etc since w0 and w1 coincide
    ws = ws(:, 2:end)';
end
