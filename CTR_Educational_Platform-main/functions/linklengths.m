% * Calculate `ll' (vector of link lengths), `kl' (vector of link curvatures) and `phil' (vector of link rotations) here *
function ll = linklengths(d, rho)
    % * YOUR CODE HERE *
    [transition_zones, ~] = sort([rho, d+rho]);
    ll = diff(transition_zones);
end