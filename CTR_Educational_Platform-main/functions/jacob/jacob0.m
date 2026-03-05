function [T, J] = jacob0(c)
    % YOUR CODE HERE
    
    n_links = size(c, 2)/3;
    J = zeros(6, 3*n_links);
    
    i = 1;
    T = eye(4);
    for j=1:3:size(c, 2)
        k = c(j);
        phi = c(j+1);
        l = c(j+2);
            
        J_link = jacobianSingleLink(k,phi,l);
        ad = adjoint(T);
        J(:, i:i+2) = ad * J_link;
        
        local_c = [k, phi, l];
        [new_T, ~] = arckinematics(local_c);
        T = T * new_T;
        i = 3 + i; 
    end 
end