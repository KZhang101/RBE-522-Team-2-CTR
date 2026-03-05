%% Ancillary functions
function [T,links] = fwkinematics(c)
    nLinks = length(c)/3;

    T = eye(4);
    links = cell(1,nLinks);

    for ii = 0 : nLinks - 1
        [Tii, link] = arckinematics([c(ii*3+1), c(ii*3+2), c(ii*3+3)]);
        links{ii+1} = applytransform(link, T);
        T = T * Tii;
    end
end