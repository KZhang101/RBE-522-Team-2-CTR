classdef Robot < handle

    properties
        % Save vector of tubes and size of vector
        tubes = []
        num_tubes = 0

        % Save link lengths, phi values, and kappa values vectors (1 x num
        % links)
        lls = []
        phi = []
        kappa = []
    end

    methods
        % Constructor. This creates an instance of the robot class 
        % Must pass in a vector of tubes
        function self = Robot(tubes)
            self.tubes = tubes;
            self.num_tubes = size(tubes, 2);
        end

        % Here we calculate the kinematics of a full CTR
        % Pass in the raw joint variables
        % Return transformation matrix that goes from home frame to end
        % effector frame
        % See functions below for each step
        function T = fkin(self, q_var)  
            % First we get the rho and theta avlues from q_var
            rho = get_rho_values(self, q_var); 
            theta = get_theta(self, q_var);

            % Next, we use rho to get the link lengths
            self.lls = get_links(self, rho);

            % Now we calculate the phi and kappa values
            [self.phi,self.kappa] = calculate_phi_and_kappa(self, theta, rho);

            % Finally we calculate the base to end effector transform
            T = calculate_transform(self, self.lls, self.phi, self.kappa);
        end

        % Get rho values from joint positions
        % Return rho (1 x i vector, i is num tubes)
        function rho = get_rho_values(self, q_var)
            % Here we extract rho values from the q_var vector

            % Initialize a vector to hold the result
            rho = zeros([1 self.num_tubes]);

            for i=2:self.num_tubes
                rho(i) = (q_var(i) - q_var(1)) * 10^-3;
            end

        end

        % Function to find the link lengths, in order
        % Returns link lengths (1 x j vector, where j is num links)
        function s = get_links(self, rho)
            % ********************************************************
            %                          TODO
            % ********************************************************
            d = [self.tubes.d];
            [transition_zones, ~] = sort([rho, d+rho]);
            transition_zones = unique(transition_zones, 'stable');
            s = diff(transition_zones);
            s = s(s ~= 0);
           
        end

        % Function to get theta values
        % Returns theta (1 x j vector where j is num links)
        function theta = get_theta(self, q_var)
            % Here we extract theta values from the q_var vector

            % Initialize a vector to hold the result
            theta = zeros([1 self.num_tubes]);

            for i=1:self.num_tubes
                theta(i) = deg2rad(q_var(i+self.num_tubes));
            end
        end

        function [chi,gamma] = linkcurvature(self, E, OD, ID, k, theta)
            % *YOUR CODE HERE*
            I = pi * (OD.^4 - ID.^4) / 64;
            chi = sum(E.*I.*k.*cos(theta)) / sum(E.*I);
            gamma = sum(E.*I.*k.*sin(theta)) / sum(E.*I);
        end


        % Function to calcualte phi values for two or three tube
        % configurations
        % Should return phi (1 x j vector, where j is num links)
        % and K (1 x j vector)
        function [phi,K] = calculate_phi_and_kappa(self, theta, rho)
            % ********************************************************
            %                          TODO
            % ********************************************************
            ll = self.lls;
            n_links = size(ll, 2);

            K   = zeros(1,n_links); % Link Curvatures [m^-1]
            phi = zeros(1,n_links); % Link Rotations [rad]
            phi_abs = zeros(1, n_links);

            d = [self.tubes.d];
            E = [self.tubes.E];
            OD = [self.tubes.od];
            ID = [self.tubes.id];
            k = [self.tubes.k];

            % sorted tube properties 
            total_tube_length = d+rho;
            length_list = [0, cumsum(ll)];

            for i = 1:n_links
                present_tubes = total_tube_length > length_list(i); % present tubes
                tube_idx = find(present_tubes == 1);
                straight_tubes = rho(tube_idx) > length_list(i); % check if the tube is straight 
                curved_tubes =  ~straight_tubes & present_tubes(tube_idx); % check if curved if not straight 
                k_tran = curved_tubes .* k(tube_idx); % only include present tubes
                
                    
                [chi,gamma] = self.linkcurvature(E(tube_idx), OD(tube_idx), ID(tube_idx), k_tran, theta(tube_idx));
                K(i) = sqrt(chi^2 + gamma^2);
                phi_abs(i) = atan2(gamma, chi);
                if i == 1
                    phi(i) = phi_abs(i);
                else
                    phi(i) = phi_abs(i) - phi_abs(i-1);
                end
                
            end 

        end

        function T = arckinematics(self, k, phi, l)
        %% YOUR CODE HERE
            if k > 0
                R = [cos(phi)*cos(k*l), -sin(phi), cos(phi)*sin(k*l);
                    sin(phi)*cos(k*l), cos(phi), sin(phi)*sin(k*l);
                    -sin(k*l), 0, cos(k*l)];
                
                p = [(cos(phi) * (1-cos(k*l))) / k;
                            (sin(phi) * (1-cos(k*l))) / k;
                            sin(k*l) / k];
            else 
                R = [cos(phi), -sin(phi), 0;
                    sin(phi), cos(phi), 0;
                    0 0 1];
                p = [0; 0; l;];
                
            end 
            T = [R p;
                0 0 0 1];

        end

        % Take in all robot dependent parameters (lls, phi, kappa) and
        % compelte the robot independent constant curvature kinamtatics
        % Returns a 4x4 transformation matrix from base frame to end
        % effector
        function T = calculate_transform(self, s, phi, K)
            % ********************************************************
            %                          TODO
            % ********************************************************
            T = eye(4);
            for ii = 1 : size(self.lls, 2)
                Tii = self.arckinematics(K(ii), phi(ii), s(ii));
                T = T * Tii;
            end

        end
        
        function [rho_des, theta_des] = inv_kin_dep(self, ind_params_des, q_current)
            % ind_params_des = desired [k1, phi1, l1, k2, phi2, l2, ...] column vector
            % q_current      = current actuation [rho1; rho2; rho3; theta1; theta2; theta3]
        
            maxIter = 500;
            tol     = 1e-4;
            alpha   = 0.4;              % step size
            lambda  = 0.01;             % damping (tune this!)
        
            q = q_current;              % start from current actuation
        
            for iter = 1:maxIter
                % Forward: compute current config from current actuation
                config_cur = forward_actuation_to_config(self, q(1:3), q(4:6));
        
                % Error in configuration space
                e = ind_params_des - config_cur;
        
                if norm(e) < tol
                    break;
                end
        
                % Numerical Jacobian using central differences
                f_handle = @(qq) forward_actuation_to_config(self, qq(1:3), qq(4:6));
                J_dep = central_diff_jacobian(f_handle, q, 1e-7);
        
                % Damped least squares
                J_dls = J_dep' / (J_dep * J_dep' + lambda^2 * eye(size(J_dep,1)));
                dq = alpha * J_dls * e;
        
                % Update actuation
                q = q + dq;
        
                % Optional: enforce physical constraints
                % q(1:3) = max(q(1:3), 0);             % rho >= 0
                % q(1:3) = sort(q(1:3));               % keep rho ordered if needed
            end
        
            if norm(e) >= tol
                warning('inv_kin_dep did not converge. Final config error: %.4f', norm(e));
            end
        
            rho_des  = q(1:3);
            theta_des = q(4:6);
        end


    %% Inverse Kinematics
    function ind_params = inv_kin_ind(self, traj, q_current)
        % traj      = 3×n matrix of desired positions [x; y; z] per column
        % q_current = current joint vector for initial guess
    
        nPts = size(traj, 2);
    
        % Get current actuation to initialize realistic configuration
        rho   = self.get_rho_values(q_current);
        theta = self.get_theta(q_current);
        lls   = self.get_links(rho);               % usually 5 segments
        [phi_cur, kappa_cur] = self.calculate_phi_and_kappa(theta, rho);
    
        % Determine number of segments/parameters from current state
        n_segs   = length(lls);
        n_params = 3 * n_segs;                     % e.g. 15 for 5 segments
    
        % Pre-allocate output: each column is one configuration (15×1 flattened)
        ind_params = zeros(n_params, nPts);
    
        % Build initial c_current as ROW vector (for fkin/jacob0 compatibility)
        c_current = zeros(1, n_params);
        for i = 1:n_segs
            idx = 3*(i-1) + 1;
            c_current(idx)   = kappa_cur(i);
            c_current(idx+1) = phi_cur(i);
            c_current(idx+2) = lls(i);
        end
    
        maxIter = 1000;
        tol     = 1e-3;
        lambda  = 0.1;
        alpha   = 0.5;
    
        for kk = 1:nPts
            xd = traj(:, kk);
    
            for iter = 1:maxIter
                % Forward kinematics
                [T, ~] = self.fkin(c_current);
                x = T(1:3, 4);
    
                e = xd - x;
                if norm(e) < tol
                    break;
                end
    
                % Jacobian
                J = jacob0(c_current);     % assumes row vector input
                Ja = jacoba(J, x);
    
                % Damped least squares
                J_dls = Ja' / (Ja*Ja' + lambda^2 * eye(3));
                deltaQ = alpha * J_dls * e;
    
                % Update
                c_current = c_current + deltaQ;
            end
    
            if norm(e) >= tol
                fprintf('Point %d did not converge. Final error: %.4f\n', kk, norm(e));
            end
    
            % Store as column
            ind_params(:, kk) = c_current(:);
        end
    end
    
    function config_vec = forward_actuation_to_config(self, rho, theta)
        lls = self.get_links(rho);          % should be 1×5 (or variable length)
    
        [phi, kappa] = self.calculate_phi_and_kappa(theta, rho);  % both 1×5
    
        n_segs = length(lls);
        if length(kappa) ~= n_segs || length(phi) ~= n_segs
            error('kappa/phi length does not match number of segments');
        end
    
        config_vec = zeros(3*n_segs, 1);
        for i = 1:n_segs
            idx = 3*(i-1) + 1;
            config_vec(idx)   = kappa(i);
            config_vec(idx+1) = phi(i);
            config_vec(idx+2) = lls(i);
        end
    end
    
    function J = central_diff_jacobian(f_handle, q_current, h)
        % f_handle    = @(q) forward_actuation_to_config(self, q(1:3), q(4:6))
        % q_current   = [rho; theta] column vector
        % h           = step size, e.g. 1e-7
    
        if nargin < 3
            h = 1e-7;
        end
    
        n = length(q_current);          % usually 6
        y0 = f_handle(q_current);       % base output
        m = length(y0);                 % e.g. 9 for 3 sections
    
        J = zeros(m, n);
    
        for k = 1:n
            q_plus = q_current;   q_plus(k)  = q_current(k) + h;
            q_minus = q_current;  q_minus(k) = q_current(k) - h;
    
            y_plus  = f_handle(q_plus);
            y_minus = f_handle(q_minus);
    
            J(:,k) = (y_plus - y_minus) / (2*h);
        end
    end
    
    %{
    % Example usage
    desired_pos = [x_des; y_des; z_des];
    
    % Step 1: Task-space → configuration space
    ind_params_des = self.inv_kin_ind(desired_pos, q_current);
    
    % Step 2: Configuration space → actuation space
    [rho_des, theta_des] = self.inv_kin_dep(ind_params_des, q_current);
    
    % Now send to motors: rho_des and theta_des
    %}
    end
end

