function J_a = jacoba(c,M)    
    % your code here
    [T, J] = jacob0(c);
    
    T = T * M;
    
    J_w = J(1:3, :); % check
    J_v = J(4:6, :);
    p = T(1:3, 4);
    
    p_skew = [0, -p(3), p(2);
          p(3), 0, -p(1);
          -p(2), p(1), 0];
        
    J_a = J_v - p_skew * J_w;
end
