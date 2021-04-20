%
% loads hypostatic 2d tensegrity
% 
function framework = load_exampleTensegrity

n = 6;
dim = 2;
h = 0;
x = [0;0;
    1;0;
    0.5;1;
    0.5;0.5;
    1/3;0 + h;
    2/3;0 + h];
a = zeros(n);
a(2,3) = 1; a(1,3) = 1;
a(1,4) = 1; a(3,4) = 1; a(2,4) = 1;
a(1,5) = 1; a(5,6) = 1; a(6,2) = 1;
a = a+a';

framework.pfix = [1,4];
framework.types = [1,1,-1,-1,-1,1,1,1];
framework.n = n;
framework.dim = dim;
framework.x = x;        % position of nodes
framework.lattice = [];   % lattice vectors (none here, since no periodic boundary conditions)
framework.ap = [];        % periodic adjacency matrix (none here since no periodic boundary conditions)
framework.a = a;        % adjacency matrix
framework.lengths = NaN;    % desired distance matrix (same size as a)


end