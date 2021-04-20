% calculate matrix of pairwise distances between particles in arbitary
% dimension

% created jan 25, 2018

function r = distances(x,dim)

n = length(x) / dim;  % number of particles

r = zeros(n);

for i=1:n
    xi = x(dim*i-(dim-1):dim*i);
    for j=i+1:n
        xj = x(dim*j-(dim-1):dim*j);
        r(i,j) = norm(xi-xj);
        r(j,i) = r(i,j);
    end
end