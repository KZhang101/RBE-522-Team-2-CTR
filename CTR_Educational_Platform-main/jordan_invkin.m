% ===================
% Quick Test Script
% ===================

clc; clear; 

% Define our tubes (ID, OD, r, l, d, E)
tube1 = Tube(3.046*10^-3, 3.3*10^-3, 1/17, 90*10^-3, 50*10^-3, 1935*10^6);
tube2 = Tube(2.386*10^-3, 2.64*10^-3, 1/22, 170*10^-3, 50*10^-3, 1935*10^6);
tube3 = Tube(1.726*10^-3, 1.98*10^-3, 1/29, 250*10^-3, 50*10^-3, 1935*10^6);

% TODO You will need to uncomment one of the below lines, depending if you
% are testing a 2-tube or 3-tube case
% tubes = [tube1, tube2];
tubes = [tube1, tube2, tube3];

robot = Robot(tubes);

% Current actuation (example values — replace with real ones)
q_current = [10, 50, 80, 45, -45, 45];   % [rho mm; theta deg] — adjust units/format

% Desired end-effector position (example)
desired_pos = [0.05; 0.02; 0.12];     % [x; y; z] in meters

% Step 1: Task → Configuration
fprintf('Running independent IK...\n');
ind_params_des = robot.inv_kin_ind(desired_pos, q_current);

fprintf('Desired configuration (first column):\n');
disp(ind_params_des(:,1));

% Step 2: Configuration → Actuation
fprintf('Running dependent IK...\n');
[rho_des, theta_des] = robot.inv_kin_dep(ind_params_des(:,1), q_current);

fprintf('Solved actuation:\n');
fprintf('  rho  = [%.4f  %.4f  %.4f] mm\n', rho_des*1000);   % if rho in meters
fprintf('  theta = [%.2f  %.2f  %.2f] deg\n', rad2deg(theta_des));

% Verification: forward check
config_check = robot.forward_actuation_to_config(rho_des, theta_des);
error_check = norm(ind_params_des(:,1) - config_check);
fprintf('Verification error (should be small): %.2e\n', error_check);




