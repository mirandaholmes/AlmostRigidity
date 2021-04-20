% Compute the stress matrix for a given stress

% Set dim=1 to send back reduced version, with no kronecker product yet

% created August 8, 2018

% updated July 9, 2019: I was computing the wrong sign! Changed sign of M0.


function M = stressmatrix(w,a,ap,dim)

n = size(a,1);
[rr,cc] = find(triu(a));
[rp,cp] = find(ap);
nb = length(rr);
nbp = length(rp);
rt = [rr;rp];
ct = [cc;cp];
nbt = nb+nbp;

M0 = sparse([rt;ct;[1:n]'], [ct;rt;[1:n]'], [w;w;zeros(n,1)],n,n);
for k=1:n
    M0(k,k) = -sum(M0(:,k));
end
M0 = -M0;   % fix july 9, 2019
% sums up stresses that correspond to same bonds, if there is a regular
% one and a periodic one, for e.g.
% i.e. quadratic form is (wk+wl)|vi-vj| if pair (i,j) is in bond k and
% bond l

% full stress matrix; expand in all particle coordinates
M = kron(M0,eye(dim));  



