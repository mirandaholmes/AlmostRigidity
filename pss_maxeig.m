% Run pss test to maximize minimum eigenvalue of stress matrix on 
% space of flexes

% created jan 29, 2018

% set ifv = 0 to find stress matrix with largest min eigenvalue



function [X,t,alph,Mevals] = pss_maxeig(V,W,dim,a,ap,ifquiet)

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

% If there are no flexes, return
if(nv == 0)
    return;
end

% If there are no stresses, return
if(nw == 0)
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


% If there is only one stress, return
if(nw == 1) 
    lam1 = min(Mevals);
    lam2 = max(Mevals);
    
    % eigenvalues are different signs, so return; no psd solution
    if( (lam1 < 0-tol && lam2 > 0+tol) || (lam2 < 0-tol && lam1 > 0+tol) )
        return;
    end
    
    % There is a psd solution; determine whether to use +M or -M
    if(lam1 > 0)
        t = lam1;
        X = M;
        alph = 1;
    else
        t = -lam2;
        X = -M;
        alph = -1;
    end
    return;
end


% If there is only one flex, return
if(nv == 1)
    [t,ind] = max(abs(M));
    X = M(ind);
    alph(ind) = 1;
    if(M(ind) < 0)
        X = -M(ind);
        alph(ind) = -1;
    end
    return;
end


% If all eigenvalues are small, return
%%% DO SOMETHING (?) ****
% **what is tolerance?**



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
A = kron(alph,eye(nv))';
variable t;
maximize t;         %  maximize minimum eigenvalue = t
y : X == M*A - t*eye(nv);
z : X == semidefinite(nv);
norm(alph) <= 1;
cvx_end

% Set X equal to linear combination of stress matrices
X = X + t*eye(nv);

% another way
%R = randn(nv)
%maximize trace(R*X);











    
    
    