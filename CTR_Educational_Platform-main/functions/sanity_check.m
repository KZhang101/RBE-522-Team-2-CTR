%% Create a robot with 2 tubes
nTubes = 2;
OD     = [1.80, 1.60] * 1e-3; % tube outer diameters [m]
ID     = [1.60, 1.40] * 1e-3; % tube inner diameters [m]
k      = [50, 100];           % tube precurvatures [m^-1]
d      = [10, 10]     * 1e-3; % length of curved section [m]
E = 75e9 * ones(1, nTubes);   % Young's modulus of Nitinol [N/m^2]

%% Set the joint variables
rho = [20, 50] * 1e-3;  % tube translations [m]
alpha = [pi/4, -pi/4];    % tube rotations [rad]

% % q values for testing with two tubes
% q_var = [0, 0, 0, 0; 
%          20, 50, 45, -45; 
%          30, 35, -35, 20];

%% Calculate the forward kinematics
% Step 1: Robot-dependent kinematics
% This robot has n = 2 tubes, therefore it will have 2*n-1 = 3 curved links
nLinks = 2*nTubes - 1;
ODlink = [OD(1) OD(1) OD(2)];

ll   = zeros(1,nLinks); % Link Lengths [m]
kl   = zeros(1,nLinks); % Link Curvatures [m^-1]
phil = zeros(1,nLinks); % Link Rotations [rad]

ll = linklengths(d, rho);
n_links = size(ll, 2);

% sorted tube properties 
total_tube_length = d+rho;
length_list = [0, cumsum(ll)];

for i = 1:n_links
    present_tubes = total_tube_length > length_list(i); % present tubes
    tube_idx = find(present_tubes == 1);
    straight_tubes = rho(tube_idx) > length_list(i); % check if the tube is straight 
    curved_tubes =  ~straight_tubes & present_tubes(tube_idx); % check if curved if not straight 
    k_tran = curved_tubes .* k(tube_idx); % only include present tubes
    
        
    [chi,gamma] = linkcurvature(E(tube_idx), OD(tube_idx), ID(tube_idx), k_tran, alpha(tube_idx));
    kl(i) = sqrt(chi^2 + gamma^2);
    if i == 1
        phil(i) = atan(gamma / chi); % delta phi 
    else
        phil(i) = atan(gamma / chi) - phil(i-1); % delta phi 
    end 
    
end 

% Step 2: Robot-independent kinematics
%% 
disp(kl)
disp(phil)
disp(ll)
T = eye(4);
links = cell(1,nLinks);

for ii = 1 : nLinks
    Tii = arckinematics(kl(ii), phil(ii), ll(ii));
    T = T * Tii;
    disp(T)
end