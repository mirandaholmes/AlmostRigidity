% compute flexes and stresses whose singular values are less than a given
% threshhold tolerance
%
% remember: 
%       # flex - # s.s. = nvar - ncon
%

% created jan 25, 2018

%
% NOTE: not getting ghost values for ifsparse = 1;
% could do QR to get ghost singular vectors
%
% NOTE 2: not getting singular values (or V) correct, if actually want to restrict
% to a particular subspace. All singular values are calculated assuming
% trivials are used to form linear space. 
%
% compute_flexstress_0: problem with singular values, as described in note
% 2 above. 
%
% modify may 6, 2019
% compute_flexstress: uses [R;L] to compute null space and ... 
% NO! problem then with stresses: get stresses associated with extra linear
% constraints. (Not if they're trivial...) 
% 
%

function [V,W,s,gap] = compute_flexstress(R,L,p)

ncon = size(R,1);  % number of constraints
nvar = size(R,2);  % number of variables

tols = p.tols;
ifsparse = p.ifsparse;
ns = p.ns;

V = [];
W = [];
gap = NaN;


% calculations in full arithmetic
if(~ifsparse)
  
    [Wf,Sf,Vf] = svd(full(R));

end

% calculations in sparse arithmetic
% *** NOTE: missing ghost singular vectors ***
if(ifsparse)
    tic
    [Wf,Sf,Vf] = svds(R,p.ns,'smallest');
    toc  
end

% singular values (decreasing order)
svals = diag(Sf);

% find index of largest singular value less than threshhold
itol = find(svals > tols,1,'last') + 1;

% find index of starting singular value, to keep only ns of them
if(~ifsparse)
    istart = max(itol, nvar-ns+1);  % index of all sing vals < tol, or ns smallest ones
else
    istart = max(itol,1); % computed only ns singular values anyways
end

% extract singular values, right sing vecs, left sing vecs
s = flipud(svals);
V1 = [];
if(~ifsparse)
    if(istart <= nvar)
        V1 = Vf(:,istart:nvar);
        V1 = fliplr(V1);
    end
    if(istart <= ncon)
        W = Wf(:,istart:ncon);
        W = fliplr(W);
    end
end
if(ifsparse && istart <= ns)
    V1 = Vf(:,istart:end);   % right null space  %V1 = Vf(:,istart:nvar); % used for Full calcs
    V1 = fliplr(V1);
    W = Wf(:,istart:end); % left null space %W = Wf(:,istart:ncon); % used for Full calcs
    W = fliplr(W);
end

% remove flexes which lie in subspace L
if(~isempty(L) && ~isempty(V1))
    A = null(L*V1);
    V = V1*A;  % solves LV2 = 0 and lives in subspace V1
    %V = orth(V);   % find an orthogonal basis of V2, via svd
    %%% unnecessary, since A is orthogonal, V1 is orthogonal
elseif(~isempty(V1))
    V = V1;
else 
    V = [];
end
%%% do something to account for singular values in L's subspace?****



% get gap
if(itol > 1 && itol <= length(svals))
    gap = svals(itol-1) - svals(itol);
else
    gap = NaN;
end



% % extract singular values, right sing vecs, left sing vecs
% s = flipud(svals); % return all singular values %s = svals(istart:nvar);
% V1 = Vf(:,istart:end);   % right null space  %V1 = Vf(:,istart:nvar); % used for Full calcs
% V1 = fliplr(V1);
% W = Wf(:,istart:end); % left null space %W = Wf(:,istart:ncon); % used for Full calcs 
% W = fliplr(W);






