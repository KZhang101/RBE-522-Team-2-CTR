% =========================
% test_ctr_ik.m
% =========================
clear; clc; close all;

% ------------------------------------------------------------
% 1) Build your tubes and robot
% ------------------------------------------------------------
% TODO: create your tube objects/structs with fields:
%   d, E, od, id, k
% and assemble into a 1xN vector "tubes"
%
% Example placeholders (REPLACE):
% tubes(1) = struct('d', ..., 'E', ..., 'od', ..., 'id', ..., 'k', ...);
% tubes(2) = struct('d', ..., 'E', ..., 'od', ..., 'id', ..., 'k', ...);
% tubes(3) = struct('d', ..., 'E', ..., 'od', ..., 'id', ..., 'k', ...);
% Define our tubes (ID, OD, r, l, d, E)
tube1 = Tube(3.046*10^-3, 3.3*10^-3, 1/17, 90*10^-3, 50*10^-3, 1935*10^6);
tube2 = Tube(2.386*10^-3, 2.64*10^-3, 1/22, 170*10^-3, 50*10^-3, 1935*10^6);
tube3 = Tube(1.726*10^-3, 1.98*10^-3, 1/29, 250*10^-3, 50*10^-3, 1935*10^6);

% TODO You will need to uncomment one of the below lines, depending if you
% are testing a 2-tube or 3-tube case
% tubes = [tube1, tube2];
tubes = [tube1, tube2, tube3];
robot = Robot(tubes);
n = robot.num_tubes;



% ------------------------------------------------------------
% 3) IK options
% ------------------------------------------------------------
opts.maxIters    = 2000;
opts.tol         = 1e-5;   % TODO: units must match tip_position() output (likely meters)
opts.lambda      = 1e-2;
opts.alpha       = 0.4;
opts.h_rho_mm    = 0.3;    % mm perturbation
opts.h_theta_deg = 0.8;    % deg perturbation
opts.q_min = [20; 20; 20;  -180; -180; -180];  % mm, deg
opts.q_max = [50; 50; 50;   180;  180;  180];  % mm, deg

% Optional bounds (strongly recommended)
% TODO: set realistic bounds for your robot
% opts.q_min = [ ... ; ... ];  % (2n)x1
% opts.q_max = [ ... ; ... ];  % (2n)x1

% ------------------------------------------------------------
% 4) Generate 5 feasible targets by random sampling with q3 > q2 > q1
% ------------------------------------------------------------
K = 5;
targets = zeros(3, K);
q_seed  = zeros(2*n, K);

% TODO: set realistic bounds for your robot/assignment
q1_min_mm = 20;      % base insertion for tube 1 (mm)
q3_max_mm = 50;      % max insertion for tube 3 (mm)

theta_min_deg = -180;
theta_max_deg = 180;

for k = 1:K
    % Sample ordered insertions: q3 > q2 > q1
    q1 = q1_min_mm + (q3_max_mm - q1_min_mm) * rand();
    q2 = q1 + (q3_max_mm - q1) * rand();   % ensures q2 >= q1
    q3 = q2 + (q3_max_mm - q2) * rand();   % ensures q3 >= q2

    rho_rand = [q1; q2; q3];

    % Sample random rotations (degrees)
    theta_rand = theta_min_deg + (theta_max_deg - theta_min_deg) * rand(n,1);

    qk = [rho_rand; theta_rand];

    xk = robot.tip_position(qk);

    targets(:,k) = xk;
    q_seed(:,k)  = qk;
end

% ------------------------------------------------------------
% 2) Initial guess q0 = [rho(mm); theta(deg)]
% ------------------------------------------------------------
rho0_mm    = zeros(n,1);   % TODO: pick feasible start insertions (mm)
theta0_deg = zeros(n,1);   % TODO: pick feasible start rotations (deg)
q0 = q_seed(:,k) + [randn(n,1)*5; randn(n,1)*10]; % 1mm, 5deg noise

% ------------------------------------------------------------
% 5) Run IK and show convergence
% ------------------------------------------------------------
infos = cell(1, size(targets,2));
q_sols = zeros(2*n, size(targets,2));

figure; hold on; grid on; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('CTR IK: Desired targets and achieved tips');

% Desired targets
scatter3(targets(1,:), targets(2,:), targets(3,:), 60, 'filled');

sum = 0;
for k = 1:size(targets,2)
    
    xd = targets(:,k);
    q0 = q_seed(:,k) + [randn(n,1)*5; randn(n,1)*10]; % 1mm, 5deg noise


    [q_sol, info] = robot.ik_nr(q0, xd, opts);

    infos{k} = info;
    q_sols(:,k) = q_sol;

    x_final = info.x_hist(:,end);
    scatter3(x_final(1), x_final(2), x_final(3), 90, 's'); % achieved tip

    fprintf('Target %d: converged=%d, iters=%d, final_err=%.6g\n', ...
        k, info.converged, info.iters, info.e_hist(end));
    sum = sum + info.e_hist(end);
end

% Calculate the average error across all targets
averageError = sum / size(targets, 2);
fprintf('Average final error across all targets: %.6g\n', averageError);

legend('Desired targets','Achieved tips');

% ------------------------------------------------------------
% 6) Plot error history for each target
% ------------------------------------------------------------
figure; hold on; grid on;
xlabel('Iteration'); ylabel('||e||');
title('IK convergence (tip position error norm)');

for k = 1:numel(infos)
    plot(infos{k}.e_hist);
end
legend('Target1','Target2','Target3','Target4','Target5');