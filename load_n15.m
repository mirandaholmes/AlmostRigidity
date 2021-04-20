
function framework = load_n15

% dimension parameters
n = 15;     % number of particles
dim = 3;    % dimension of framework

% coordinates
x = [
0.0000, 0.0000, 0.0000,...
0.0022, 1.4459, 1.3818,...
0.0022, -1.4459, -1.3818,...
1.3825, 0.0000,  1.4452,...
1.3825, 0.0000, - 1.4452,...
1.4466,  1.3810,  0.0000,...
1.4466, -1.3810, - 0.0000,...
-0.5235, -0.4171, 1.8847,...
-0.5235, 0.4171,-1.8847,...
-0.4166, 1.8846, -0.5240,...
-0.4166, -1.8846, 0.5240,...
- 1.7301, 0.6859, -0.7324,...
- 1.7301, -0.6859, 0.7324,...
-2.1985, -1.2156, 1.1384,...
-2.1985, 1.2156, -1.1384]';

% adjacency matrix 
tolD = 2e-4;  % tolerance for saying particles are bonded
a = get_adj(x,2+tolD,dim);


% set up framework struct
framework.x = x;        % position of nodes
framework.pfix = [];      % particles to fix (none here; use COM constraints)
framework.lattice = [];   % lattice vectors (none here, since no periodic boundary conditions)
framework.ap = [];        % periodic adjacency matrix (none here since no periodic boundary conditions)
framework.a = a;        % adjacency matrix
framework.lengths = 2*a;    % desired distance matrix (same size as a)
framework.n = n;        % # of nodes
framework.dim = dim;      % dimension framework 
framework.types = [];     % types, if tensegrity (0=bar, 1=strut, -1=cable)
                          % empty means all are bars
                          
end
