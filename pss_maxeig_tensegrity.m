% Run pss test to check if tensegrity is prestress stable
%
% types: 1 = cable, 0 = bar, -1 = strut. 
% careful -- order in types, must match order of bonds in [rr,cc] = find(triu(a)). 
%
%
% created july 8, 2019 from pss_maxeig.m
%


function [X,t,alph,Mevals] = pss_maxeig_tensegrity(V,W,types,dim,a,ap,ifquiet)

tol = 8*eps;  % tolerance for 0


nv = size(V,2);
nw = size(W,2);

% default values to return
X = [];
Mevals = [];
alph = zeros(1,nw);
t = NaN;
y = NaN;
z = NaN;


% If there are no stresses, return
if(nw == 0)
    return;
end
 
% Construct E-matrix
E = diag(types);
E(types == 0,:) = []; % it only acts on cables or struts
nt = size(E,1);


% If there are no flexes, look for stress that blocks all tensegrity elements

if(nv == 0)
    if(ifquiet)
        cvx_begin quiet
    else
        cvx_begin
    end
    variable alph(1,nw)
    variable s
    dual variable u
    maximize s;         
    u : E * W*alph' - s*ones(nt,1) >= 0;
    norm(alph) <= 1;
    cvx_end
    
    if(norm(alph) < tol || s < tol)
        alph = NaN;
    end
    return;
end


% Construct stress matrices
M = zeros(nv,nv*nw);   % holds all stress matrices
Mevals = zeros(nv,nw);  % keep track of eigenvalues of stress matrices
for jw=1:nw
    
    M1 = stressmatrix(W(:,jw),a,ap,dim);
    
    % work in V-space
    Mt = V'*M1*V;
    Mt = (Mt + Mt')/2;  % symmetrize, to fix numerical issues with imaginary eigenvalues
    M(:,(jw-1)*nv+1:jw*nv) = Mt;
    
    % record eigenvalues
    Mevals(:,jw) = eig(Mt);
end


% Run cvx to find the stress matrix with the maximum minimum eigenvalue. 

if(ifquiet)
    cvx_begin quiet 
else 
    cvx_begin
end
variable alph(1,nw)
variable X(nv,nv) 
dual variable y
dual variable z
dual variable u
A = kron(alph,eye(nv))';
variable t;
variable s;
maximize min(t,s);         %  maximize minimum eigenvalue = t
y : X == M*A - t*eye(nv);
z : X == semidefinite(nv);
u : E * W * alph'  - s*ones(nt,1) >= 0;
norm(alph) <= 1;
cvx_end


% Set X equal to linear combination of stress matrices
X = X + t*eye(nv);



% another way
%R = randn(nv)
%maximize trace(R*X);











    
    
    