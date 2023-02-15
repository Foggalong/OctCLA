% This file is part of OctCLA, Copyright (c) 2022 Josh Fogg, released
% under the MIT License. See: https://github.com/Foggalong/OctCLA

function ws = calculate_turningpoints(mu, covar, lb, ub, KKT)
    % CALCULATE_TURNINGPOINTS return portfolios for CLA turning points
    %
    % Takes a vector of expected returns (mu), a covariance matrix
    % (covar), a vector of lower bounds on assest weights (lb), and a
    % vector of upper bounds on asset weights (ub), then calculates
    % the turning points of the corresponding portfolio using the
    % Critical Line Algorithm. The function returns a matrix of the
    % portfolio asset weights at each of those points (ws).
    %
    % Takes an optional integer argument (KKT) which indicates whether
    % calculated turning points should be verified through a KKT
    % conditions check and if so by what type of check (0: none,
    % 1: simple (default), 2: full, 3: both).
    %
    % See also, STARTING_SOLUTION, BECOMES_FREE, MOVE_TO_BOUND

    % set default value for KKT
    if (nargin < 5); KKT = 1; end

    % calculate starting solution
    [F, B, ws] = starting_solution(mu, lb, ub);

    % initial inversion, the only one calculated without shortcuts
    % TODO this is 1x1, should be able to just take reciprocal 
    invcovarF = inv(covar(F,F));

    lam_current = inf;
    t = 1;  % so current w will be ws(:,t)

    while true
        % case a where a free asset moves to its bound
        [i_ins, lam_ins, b, d_ins] = move_to_bound(mu, covar, invcovarF, lb, ub, F, B, lam_current, ws(:,t), KKT);
        % case b where an asset on its bound becomes free
        [i_outs, lam_outs, d_outs] = becomes_free(mu, covar, invcovarF, lb, ub, F, B, lam_current, ws(:,t), KKT);

        if (i_ins ~= NaN || i_outs ~= NaN)
            lam_current = max(lam_ins, lam_outs);
            % if lam < 0 the risk is increasing again
            % make lam = 0 the last iteration
            if lam_current < 0; lam_current = 0; end

            % add an additional column to ws for the new asset weights
            ws = [ws, zeros(length(mu),1)];
            ws(B, t+1) = ws(B, t);
            t = t+1;

            % update gamma  % TODO make this neater
            if isempty(B)
                gam_top = -lam_current*sum(invcovarF)*mu(F) + 1;
            else
                gam_top = -lam_current*sum(invcovarF)*mu(F) + 1 - sum(ws(B,t)) + sum(invcovarF)*(covar(F,B)*ws(B,t));
            end
            gam = gam_top/sum(sum(invcovarF));

            % update w_F^(t)  % TODO make this neater
            if isempty(B)
                ws(F,t) = gam*sum(invcovarF,2) + lam_current*(invcovarF*mu(F));
            else
                ws(F,t) = -invcovarF*(covar(F,B)*ws(B,t)) + gam*sum(invcovarF,2) + lam_current*(invcovarF*mu(F));
            end

            % if lambda = 0 then risk is increasing again; terminate
            if lam_current <= 0; break; end

            % update free and bounded asset index sets
            if (lam_ins > lam_outs)
                % update the inverse
                % need index in inverse, not full matrix
                j = find(F==i_ins);  % TODO could find be replaced?
                invcovarF = inverse_shrink(invcovarF, j);
                % bound weight i_ins
                F = F(F ~= i_ins);  % F = F\{i}
                B = [B, i_ins];     % B = Bu{i}
                ws(i_ins, t) = b;   % w_i_inside^(t) = b  % TODO is this needed?
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
            % NOTE doesn't check lam=0 since not (necessarily) a TP
            if (KKT > 0)
                disp('Checking KKT conditions...')
                errors = 0;

                if (KKT == 1) || (KKT == 3) 
                    errors = errors + ~kkt_full(lb, ub, d, ws(:,t), covar, lam_current, mu, gam);
                end

                if (KKT == 2) || (KKT == 3)
                    errors = errors + ~kkt_partitioned(covar, ws(:,t), lam_current, gam, mu, B, F);
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
