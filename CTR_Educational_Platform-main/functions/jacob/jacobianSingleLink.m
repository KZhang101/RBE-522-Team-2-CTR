function J = jacobianSingleLink(k,phi,l)
    % YOUR CODE HERE
    J = zeros(6, 3);
    J(4, 1) = cos(phi) * (cos(k * l) - 1) / k^2;
    J(5, 1) = sin(phi) * (cos(k * l) - 1) / k^2;
    J(6, 1) = -(sin(k * l) - k * l) / k^2;
    J(1, 1) = -l * sin(phi);
    J(2, 1) = l * cos(phi);
    
    J(3, 2) = 1;
    
    J(6, 3) = 1;
    J(1, 3) = -k * sin(phi);
    J(2, 3) = k * cos(phi);
end