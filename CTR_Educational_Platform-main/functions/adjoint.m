function Vtrans = adjoint(T)
        % your code here
        R = T(1:3, 1:3);
        p = T(1:3, 4);
        adj = eye(6);
        adj(1:3, 1:3) = R; 
        adj(4:6, 4:6) = R;
        adj(4:6, 1:3) = skew(p') * R;
        Vtrans = adj;
end