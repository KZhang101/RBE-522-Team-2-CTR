% =========================
% Robot.m
% =========================
classdef Robot < handle

    properties
        % Save vector of tubes and size of vector
        tubes = []
        num_tubes = 0

        % Save link lengths, phi values, and kappa values vectors (1 x num links)
        lls   = []
        phi   = []
        kappa = []
    end

    methods
        % Constructor. This creates an instance of the robot class.
        % Must pass in a vector (1 x num_tubes) of tube objects/structs
        function self = Robot(tubes)
            self.tubes = tubes;
            self.num_tubes = size(tubes, 2);
        end

        % ==========================================================
        % Forward Kinematics (motor space -> tip pose)
        % q_var convention (as you already use):
        %   q_var = [rho_1 ... rho_n, theta_1 ... theta_n]^T
        %   rho_i in mm  (relative insertions are computed in get_rho_values)
        %   theta_i in degrees
        % Returns 4x4 transform base->tip
        % ==========================================================
        function T = fkin(self, q_var)
            % 1) Extract rho and theta from q_var
            rho   = self.get_rho_values(q_var);
            theta = self.get_theta(q_var);

            % 2) Link lengths from rho
            self.lls = self.get_links(rho);

            % 3) Compute phi and kappa along links
            [self.phi, self.kappa] = self.calculate_phi_and_kappa(theta, rho);

            % 4) Build transform using constant curvature arc kinematics
            T = self.calculate_transform(self.lls, self.phi, self.kappa);
        end

        % Tip position only (3x1)
        function x = tip_position(self, q_var)
            T = self.fkin(q_var);
            x = T(1:3,4);
        end

        % ==========================================================
        % Helpers: pack/unpack (optional convenience)
        % ==========================================================
        function q = pack_q(self, rho_mm, theta_deg)
            q = [rho_mm(:); theta_deg(:)];
        end

        function [rho_mm, theta_deg] = unpack_q(self, q)
            n = self.num_tubes;
            rho_mm    = q(1:n);
            theta_deg = q(n+1:2*n);
        end

        % ==========================================================
        % Get rho values from joint positions
        % Return rho (1 x num_tubes)
        %
        % NOTE: This function encodes your modeling choice:
        % rho(i) = (q_var(i) - q_var(1))*1e-3, i=2..n
        % rho(1) stays 0 (base reference)
        % ==========================================================
        function rho = get_rho_values(self, q_var)
            rho = zeros(1, self.num_tubes);
            for i = 1:self.num_tubes
                rho(i) = q_var(i) * 1e-3; % mm -> m (absolute insertion from base)
            end
        end

        % ==========================================================
        % Function to find the link lengths, in order
        % Returns link lengths (1 x num_links)
        % ==========================================================
        function s = get_links(self, rho)
            d = [self.tubes.d]; % curved lengths (m)
        
            % Transition zones include base 0
            transition_zones = [0, rho, rho + d];
        
            % Keep only nonnegative zones (base reference)
            transition_zones = transition_zones(transition_zones >= 0);
        
            transition_zones = unique(sort(transition_zones), 'stable');
        
            s = diff(transition_zones);
        
            % Remove tiny/zero segments
            s = s(s > 1e-9);
        
            % Safety: if everything got filtered, force a tiny straight segment
            if isempty(s)
                s = 1e-6;
            end
        end

        % ==========================================================
        % Function to get theta values
        % Returns theta (1 x num_tubes) in radians
        % ==========================================================
        function theta = get_theta(self, q_var)
            theta = zeros([1 self.num_tubes]);
            for i = 1:self.num_tubes
                theta(i) = deg2rad(q_var(i + self.num_tubes));
            end
        end

        % ==========================================================
        % Link curvature composition helper
        % ==========================================================
        function [chi, gamma] = linkcurvature(self, E, OD, ID, k, theta)
            % Safe composite curvature calculation
            I = pi * (OD.^4 - ID.^4) / 64;
        
            denom = sum(E .* I);
            if isempty(denom) || denom < 1e-12 || ~isfinite(denom)
                chi = 0;
                gamma = 0;
                return;
            end
        
            chi   = sum(E .* I .* k .* cos(theta)) / denom;
            gamma = sum(E .* I .* k .* sin(theta)) / denom;
        
            if ~isfinite(chi), chi = 0; end
            if ~isfinite(gamma), gamma = 0; end
        end

        % ==========================================================
        % Calculate phi (link-relative) and K (link curvature magnitudes)
        % ==========================================================
        function [phi, K] = calculate_phi_and_kappa(self, theta, rho)
            ll = self.lls;
            n_links = size(ll, 2);

            K       = zeros(1, n_links); % [m^-1]
            phi     = zeros(1, n_links); % [rad] relative per-link
            phi_abs = zeros(1, n_links); % [rad] absolute per-link

            d  = [self.tubes.d];
            E  = [self.tubes.E];
            OD = [self.tubes.od];
            ID = [self.tubes.id];
            k  = [self.tubes.k];

            total_tube_length = d + rho;      % per-tube end position
            length_list = [0, cumsum(ll)];    % link boundary positions

            for i = 1:n_links
                present_tubes = total_tube_length > length_list(i);
                tube_idx = find(present_tubes == 1);

                if isempty(tube_idx)
                    K(i) = 0;
                    if i == 1
                        phi_abs(i) = 0;
                        phi(i) = 0;
                    else
                        phi_abs(i) = phi_abs(i-1);
                        phi(i) = 0;
                    end
                    continue;
                end

                s0 = length_list(i);                         % link start (m)
                straight_tubes = (s0 < rho(tube_idx));       % before curved section begins
                curved_tubes   = ~straight_tubes;            % in curved section for those present

                k_tran = curved_tubes .* k(tube_idx); % only curved + present tubes contribute curvature

                [chi, gamma] = self.linkcurvature( ...
                    E(tube_idx), OD(tube_idx), ID(tube_idx), k_tran, theta(tube_idx));

                K(i)       = sqrt(chi^2 + gamma^2);
                phi_abs(i) = atan2(gamma, chi);

                if i == 1
                    phi(i) = phi_abs(i);
                else
                    phi(i) = phi_abs(i) - phi_abs(i-1);
                end
            end
        end

        % ==========================================================
        % Constant-curvature arc kinematics for one link
        % ==========================================================
        function T = arckinematics(self, k, phi, l)
            epsk = 1e-6;
            if abs(k) > epsk
                R = [ cos(phi)*cos(k*l), -sin(phi),  cos(phi)*sin(k*l);
                      sin(phi)*cos(k*l),  cos(phi),  sin(phi)*sin(k*l);
                      -sin(k*l),          0,         cos(k*l)];

                p = [ (cos(phi) * (1 - cos(k*l))) / k;
                      (sin(phi) * (1 - cos(k*l))) / k;
                      sin(k*l) / k ];
            else
                R = [ cos(phi), -sin(phi), 0;
                      sin(phi),  cos(phi), 0;
                      0,         0,        1];
                p = [0; 0; l];
            end

            T = [R p;
                 0 0 0 1];
        end

        % ==========================================================
        % Compose full transform across all links
        % ==========================================================
        function T = calculate_transform(self, s, phi, K)
            T = eye(4);
            for ii = 1:size(s, 2)
                Tii = self.arckinematics(K(ii), phi(ii), s(ii));
                T = T * Tii;
            end
        end

        % ==========================================================
        % CENTRAL-DIFFERENCE Jacobian in MOTOR SPACE:
        % J = d x_tip / d q_var
        %
        % q_var uses your convention:
        %   first n entries: rho in mm
        %   next n entries: theta in degrees
        %
        % h_rho_mm, h_theta_deg are small perturbations in those units.
        % ==========================================================
        function J = jacobian_cd(self, q_var, h_rho_mm, h_theta_deg)
            n = self.num_tubes;
            J = zeros(3, 2*n);

            for i = 1:(2*n)
                dq = zeros(size(q_var));

                if i <= n
                    dq(i) = h_rho_mm;      % mm
                else
                    dq(i) = h_theta_deg;   % deg
                end

                x_plus  = self.tip_position(q_var + dq);
                x_minus = self.tip_position(q_var - dq);

                if any(~isfinite(x_plus)) || any(~isfinite(x_minus))
                    J(:, i) = 0;
                else
                    J(:, i) = (x_plus - x_minus) / (2 * dq(i));
                end            
            end
        end

        % ==========================================================
        % Inverse Kinematics via Newton / Damped Least Squares
        %
        % Solves for q_var to reach desired tip position x_des (3x1).
        %
        % opts fields you should set in your test script:
        %   opts.maxIters
        %   opts.tol           (same units as x_des/x_tip; likely meters)
        %   opts.lambda
        %   opts.alpha
        %   opts.h_rho_mm
        %   opts.h_theta_deg
        % Optional:
        %   opts.q_min, opts.q_max  (2n x 1 bounds)
        %
        % Returns:
        %   q_sol
        %   info.q_hist, info.x_hist, info.e_hist, info.converged
        % ==========================================================
        function [q_sol, info] = ik_nr(self, q_init, x_des, opts)

            q = q_init(:);
            x_des = x_des(:);

            maxIters = opts.maxIters;
            tol      = opts.tol;
            lambda   = opts.lambda;
            alpha    = opts.alpha;

            q_hist = zeros(numel(q), maxIters);
            x_hist = zeros(3, maxIters);
            e_hist = zeros(1, maxIters);

            for it = 1:maxIters
                x = self.tip_position(q);
                e = x_des - x;

                q_hist(:, it) = q;
                x_hist(:, it) = x;
                e_hist(it)    = norm(e);

                if norm(e) < tol
                    break;
                end

                J = self.jacobian_cd(q, opts.h_rho_mm, opts.h_theta_deg); % 3 x (2n)

                % Damped least squares (stable Newton step)
                A  = J*J' + (lambda^2)*eye(3);

                y = A \ e;
                if any(~isfinite(y))
                    y = pinv(A) * e;
                end
                dq = J' * y;


                q = q + alpha * dq;

                % Optional bounds
                if isfield(opts, 'q_min') && isfield(opts, 'q_max')
                    q = max(q, opts.q_min(:));
                    q = min(q, opts.q_max(:));
                end
            end

            k = it;
            q_sol = q;

            info.q_hist = q_hist(:, 1:k);
            info.x_hist = x_hist(:, 1:k);
            info.e_hist = e_hist(1:k);
            info.iters  = k;
            info.converged = (e_hist(k) < tol);
        end

    end
end