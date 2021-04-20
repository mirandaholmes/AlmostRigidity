% Compute kappa given R, stress matrix, lambda

% created August 8, 2018

function kappa = pss_kappa(R,M,L,lam,ifquiet)

sr= size(R,2);
sx = sr - size(L,1);

% Compute action of R'*R, M on space orthogonal to trivial motions
N = null(full(L));  % orthogonal basis of space perpendicular to L
RRt = N'*(R'*R)*N;
Mt = N'*M*N;

% Run cvx to find the stress matrix with the maximum minimum eigenvalue. 
if(ifquiet)
    cvx_begin quiet 
else 
    cvx_begin
end
variable kappa
variable X(sx,sx) 
dual variable y
dual variable z
minimize kappa;
y : X == kappa*2*RRt + Mt - lam*eye(sx);
z : X == semidefinite(sx);

cvx_end