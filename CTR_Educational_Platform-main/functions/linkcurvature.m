function [chi,gamma] = linkcurvature(E, OD, ID, k, theta)
    % *YOUR CODE HERE*
    I = pi * (OD.^4 - ID.^4) / 64;
    chi = sum(E.*I.*k.*cos(theta)) / sum(E.*I);
    gamma = sum(E.*I.*k.*sin(theta)) / sum(E.*I);
end