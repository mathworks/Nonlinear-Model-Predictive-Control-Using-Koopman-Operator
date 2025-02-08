% Create the MyNMPC class implementation
classdef MyNMPC
    methods (Static)
        function x_next = StateFcn(x, u, pvstate)
            % x is now [ x1; x2; x3 ], with x3 = integral of error
            % u is the control input
            % pvstate(1) = reference for x1

            persistent obsFcn K

            if isempty(K)
                koopman = coder.load("KoopmanMatrix.mat", "K");
                K = koopman.K;
                obsFcn = @(x1,x2,u,r) {
                    x1;                  % 1 State 1
                    x2;                  % 2 State 2
                    r;                   % 3 Reference
                    r - x1;              % 4 Tracking error
                    x2^2;                % 5 Kinetic energy-like
                    u;                   % 6 Control
                    u^2;                 % 7 Control energy
                    (r - x1)^2;          % 8 Squared tracking error
                    x1*x2;               % 9 Mixed term
                    (r - x1)*x2;         % 10 Error-velocity product
                    };
            end

            % Unpack the augmented state
            x1 = x(1);
            x2 = x(2);
            xi = x(3);  % integral of error

            % Reference
            ref = pvstate(1);

            % ============== 1) Koopman update for [x1, x2] ==============
            currentObs = obsFcn(x1, x2, u, ref);
            psi_current = cell2mat(currentObs);

            psi_next = K * psi_current;

            x1_next = psi_next(1);
            x2_next = psi_next(2);

            % ============== 2) Integral-of-error update ==============
            % xI(k+1) = xI(k) + [ref - x1(k)]
            % (You can multiply by dt if you want "true" integral)
            xi_next = xi + (ref - x1);

            % ============== Combine for next state vector ==============
            x_next = [ x1_next; x2_next; xi_next ];
        end

        function cost = CostFcn(stageNumber, x, u, pvstage)
            % x is [ x1; x2; xI ]
            % pvstage = [ref; predHorizon; q; r; s; qi]                     

            ref = pvstage(1);
            predHorizon = pvstage(2);
            Q = pvstage(3);
            R = pvstage(4);
            S = pvstage(5);
            Qi = pvstage(6);

            x1 = x(1);
            xi = x(3);
            e = (ref - x1);    % tracking error

            if stageNumber <= predHorizon
                % Stage cost: penalize error^2, integral^2, control^2
                cost = Q*e^2 + Qi*xi^2 + R*(u^2);
            else
                % Terminal cost: can also penalize xi or not
                cost = S*e^2 + Qi*xi^2;
            end
        end


        function [DcostDx, DcostDu] = CostJacFcn(stageNumber, x, u, pvstage)

            % Extract parameters
            ref = pvstage(1);
            predHorizon = pvstage(2);
            Q = pvstage(3);
            R = pvstage(4);  
            S = pvstage(5);
            Qi = pvstage(6);

            % Extract states
            x1 = x(1);
            xi = x(3);          % Integral of error

            e = (ref - x1);     % Tracking error

            if stageNumber <= predHorizon
                % Stage cost derivatives
                % Partial derivative w.r.t x1: -2*q*e
                dCost_dx1 = -2 * Q * e;

                % Partial derivative w.r.t x2: 0 (no dependence)
                dCost_dx2 = 0;

                % Partial derivative w.r.t xI: 2*qI*xI
                dCost_dxI = 2 * Qi * xi;

                % Combine into a [3 x 1] vector
                DcostDx = [dCost_dx1; dCost_dx2; dCost_dxI];

                % Partial derivative w.r.t u: 2*r_weight*u
                DcostDu = 2 * R * u;
            else
                % Terminal cost derivatives
                % Partial derivative w.r.t x1: -2*s*e
                dCost_dx1 = -2 * S * e;

                % Partial derivative w.r.t x2: 0 (no dependence)
                dCost_dx2 = 0;

                % Partial derivative w.r.t xI: 2*qI*xI
                dCost_dxI = 2 * Qi * xi;

                % Combine into a [3 x 1] vector
                DcostDx = [dCost_dx1; dCost_dx2; dCost_dxI];

                % Partial derivative w.r.t u: 0 (no control in terminal cost)
                DcostDu = zeros(size(u));
            end
        end

    end
end