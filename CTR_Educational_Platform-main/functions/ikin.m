function q = ikin(c, target_position)
    % currentQ: 1xn vector of initial joint variables 
    % targetPose: 6x1 twist representing the target pose
   
    maxCycles = 700;
    
    
    [T, ~] = fwkinematics(c);
    
    
    iterations = 0;
    current_c = c;

    current_position = T(1:3, end);

    while norm(current_position - target_position) > 1e-6
        J = jacoba(current_c, eye(4)); % Change here!!!!!!!!!!
        J = J(1:3, :);

        error = target_position - current_position;
    
        % Least Squares method
        lambda = .1;
        J_star = J' * pinv( (J*J' + lambda^2 * eye(3)) );
        
        % Apply alpha to change speed rate
        alpha = .1;
        deltaQ = alpha * J_star * error;
        
        current_c = current_c + deltaQ';

        [T, ~] = fwkinematics(current_c);
        current_position = T(1:3, end);
        iterations = iterations + 1; 
        norm(current_position - target_position);
        if iterations > maxCycles
            fprintf("Failed to converge\n");
            norm(current_position - target_position)
            break;
        end 
    end 

    % qlims = deg2rad([-180, 180;
    %      -125, 125;
    %      -138, 138;
    %      -270, 270;
    %      -120, 133.5;
    %      -270, 270]);

    q = current_c;

    if iterations < maxCycles
        % if all(q >= qlims(:,1) & q <= qlims(:,2))
        fprintf("pass\n")
        % else
        %     fprintf("failed q lims\n")
        %     rad2deg(q);
        %     all(q >= qlims(:,1) & q <= qlims(:,2));
        %     q = zeros(6,1);
        % end 
    else
        q = zeros(size(c, 2),1);
    end 
end 

