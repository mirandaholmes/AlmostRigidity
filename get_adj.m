% Get adjacency matrix from list of particle positions

% Created Sept 7, 2013

function a = get_adj(x,bondlength,dim)

n = length(x) / dim;
r = distances(x,dim);
a = +( (r + 2*bondlength*eye(n))< bondlength);