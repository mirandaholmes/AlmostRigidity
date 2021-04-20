%
% loads hypostatic framework of n=10 unit spheres
% 
function framework = load_n10

% dimension parameters
n = 10;     % number of particles
dim = 3;    % dimension of framework

% coordinates
x = [
                   0
                   0
                   0
   1.000000000000000
   0.000000000000000
   0.000000000000000
  -0.500000000000000
   0.866025403784439
   0.000000007289415
   1.000000003306545
   1.603750749657996
   0.453609204877056
   0.999999994048218
   0.577350265753363
  -0.816496583357531
  -0.000000003967855
   1.539600726278313
  -0.544331039863996
   0.000000003306546
   1.603750740716201
   0.453609226708918
   0.999999996032145
   1.539600715548160
  -0.544331060431297
   1.500000000000000
   0.866025403784439
  -0.000000007289415
   0.500000000000000
   0.866025403784439
  -0.000000000000000
];


% adjacency matrix 
tolD = 1e-5;  % tolerance for saying particles are bonded
a = get_adj(x,1+tolD,dim);
pfix = get_pfix(a);  % fix 3 particles  
%pfix = [];  % fix COM & rotation constraints


% set up framework struct
framework.x = x;        % position of nodes
framework.pfix = pfix;      % particles to fix (none here; use COM constraints)
framework.lattice = [];   % lattice vectors (none here, since no periodic boundary conditions)
framework.ap = [];        % periodic adjacency matrix (none here since no periodic boundary conditions)
framework.a = a;        % adjacency matrix
framework.lengths = a;    % desired distance matrix (same size as a)
framework.n = n;        % # of nodes
framework.dim = dim;      % dimension framework 
framework.types = [];     % types, if tensegrity (0=bar, 1=strut, -1=cable)
                          % empty means all are bars
                          
end