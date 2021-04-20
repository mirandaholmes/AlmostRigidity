% Form rigidity matrix R corresponding to given bond distances
% Works in arbitrary dimension
%
% Currenty doesn't accept spring constants, but could change this 
% **TO DO **
%
% Created jan 25, 2018
%

function R = formconstraints_bonds(x,a,dim,varargin)


% extract periods, if needed
if(nargin > 3)
    ap = varargin{1};
    lattice = varargin{2};
else
    ap = [];
    lattice = [];
end



[rr,cc] = find(triu(a));  % regular constraints
[rp,cp,sp] = find(ap);    % periodic constraints
nb = length(rr);
nbp = length(rp);
n = size(a,1);  % number of particles

R = sparse(nb+nbp,dim*n);  % rigidity matrix


for ji=1:nb   % bond constraints
    ir = rr(ji);
    ic = cc(ji);
    p1 = x(dim*ir-(dim-1):dim*ir);  % first point
    p2 = x(dim*ic-(dim-1):dim*ic);  % second point
    R(ji, dim*ir-(dim-1): dim*ir) = (p1 - p2)';
    R(ji, dim*ic-(dim-1): dim*ic) = -(p1 - p2)';
end

% periodic constraints (only works in 2d right now)
for ji=1:nbp
    ir = rp(ji);
    ic = cp(ji);
    p1 = x(dim*ir-(dim-1):dim*ir);  % first point
    p2 = x(dim*ic-(dim-1):dim*ic);  % second point
    
    if(sp(ji) < 0)
        p1 = p1 + lattice(:,1);
    else
        p1 = p1 + lattice(:,2);
    end
    
    R(nb+ji, dim*ir-(dim-1): dim*ir) = (p1 - p2)';
    R(nb+ji, dim*ic-(dim-1): dim*ic) = -(p1 - p2)';
end




