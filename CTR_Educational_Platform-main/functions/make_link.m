function rotLink = make_link(q)
% MAKELINK Takes as input a vector q = [k, phi, l] and returns a list of
% points representing the shape of the continuum link in space
    
    % copy the configuration
    k   = q(1);
    phi = q(2);
    l   = q(3);
   
    if k == 0
        link = zeros(3, 21);
        link(3,:) = 0 : l / 20 : l;
        rotLink = link;
        
    else
        radius = 1 / k;
        theta = 0 : l * abs(k) / 20 : l * abs(k);
        
        
        % if negative curvature, bend the other way
        if k < 0, theta = -theta; end
        
        % generate the link (homogeneous coordinates)
        link = radius .* [(1 - cos(theta));
            zeros(1, length(theta));
            sin(theta);
            ones(1, length(theta)) / radius];
        
        
        % rotate the link about Z by the angle phi
        xRot = [0 -1 0 0;
            1 0 0 0;
            0 0 0 0;
            0 0 0 0];
        
        T = expm(xRot * phi);
        
        rotLink = cellfun(@(p) T * p, num2cell(link, 1), ...
            'UniformOutput', false);
        
        rotLink = cell2mat(rotLink);
        rotLink = rotLink(1:3,:);
    end
end