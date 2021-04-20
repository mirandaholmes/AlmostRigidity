%
% Test for almost-rigidity of a framework. 
% Follows the test in 
%     Holmes‐Cerfon, Miranda, Louis Theran, and Steven J. Gortler. 
%     "Almost‐Rigidity of Frameworks." Communications on Pure and Applied Mathematics (2019).
%
% Needs CVX to be pre-installed. You can download CVX from this site:
%     http://cvxr.com/cvx/
%
%
% Created April 20, 2021.
%
%

clear;

%%%===========================%%%
%%%    Parameters & Data      %%%
%%%===========================%%%

% Load data into framework object
framework = load_n10;     % n=10 hypostatic framework
%framework = load_exampleTensegrity;  % example of a 2d tensegrity
%framework = load_n15;      % 


% Parameters for almost-rigidity test
% Depending on which framework is loaded, you may want to change 'tols'. 
params = struct(...
    'tols',1e-5,...   % tolerance on singular values for null space
    'ns',inf,...      % max # of sing vectors to compute
    'ifsparse',0,...  % use sparse arithmetic to compute singular values IMPORTANT to set =0
    'lamfac',0.5...   % fraction of max eigenvalue lam0 of PSD stress matrix (lambda = lamfac*lam0)
    );

% Parameters for displaying / plotting information in this script
ifdisp = 1;  % display information as we go along
ifquiet = 1;     % use cvx in quiet mode
ifplot = 1;   % plot framework at the end
   whichWplot = 1;  % which stress to plot: 1=W, 2=W*alph
   ivplot = inf;  % which flex to plot (increasing order of sing vals), inf for all
   iwplot = 1;    % which stress to plot  (increasing order of sing vals)
   davg = 1;     % typical distance between points (for plotting)


%%%================================%%%
%%%       Set up framework         %%%
%%%================================%%%

% extract variables for ease of reference
x = framework.x;
a = framework.a;
n = framework.n;
dim = framework.dim;
ap = framework.ap;
lattice = framework.lattice;
pfix = framework.pfix;
types = framework.types;

% center cluster
x = center(x,dim); 

% Calculate derived quantities
z = full(max(sum(a)));  % max adjacency
if(isempty(types) || isempty(find(abs(types)==1,1)))
    isbar = 1; istens = 0;   % entirely a bar framework
else
    isbar = 0; istens = 1;   % a tensegrity
end

%%%================================%%%
%%%   Form constraint matrices     %%%
%%%================================%%%
% Rigidity matrix (bonds only)
R = formconstraints_bonds(x,a,dim,ap,lattice);

nb = size(R,1);     % number of bonds
nvar = size(R,2);   % number of variables

% Additional, linear constraints: either 
% (i) fix center of mass and no rotations, or (ii) fix 2 or 3 particles
if(isempty(pfix))  % (i)
    C = [];
    Lrot = formconstraints_rot(x,dim);   % infinitesimal rotation constraints
    C = [C;Lrot];
    Lcom = formconstraints_com(x,dim);   % center of mass constraints
    C = [C;Lcom];
else     % (ii)
    C = formconstraints_pfix(x,pfix,dim);  % fix 2 or 3 particles
end



%%%=================================%%%
%%%   Compute Flexes and stresses   %%%
%%%=================================%%%

% Compute flexes and stresses and gap (based on tolerances)
% Note that compute_flexstress computes singular values of R only,
% then projects singular vectors to orthogonal space to L.
% Hence, it's not the same as computing the singular values of [R;L]
% all together. (Therefore, it might not matter if we use pfix or not,
% since we end up with the same result in the end.)
[V,W,s,gap] = compute_flexstress(R,C,params);
% Output: 
%   V = almost-flexes (col#1 = smallest singular value vector, col# = second
%                      smallest, etc)
%   W = almost-stresses (same format)
%   s = list of params.ns smallest singular values 
%   gap = difference between largest singular value less than params.tols,
%         and smallest singular value greater than params.tols
%

nv = size(V,2);   % number of almost-flexes
nw = size(W,2);   % number of almost-stresses

% Display diagnostic information
if(ifdisp)
    disp('---  Spectral properties of R:  ---');
    display(['nbond = ',num2str(nb),', nvar = ',num2str(nvar)]);
    display(['nv = ',num2str(nv),', nw = ',num2str(nw)]);
    display(['zmax  = ',num2str(z)]);
    display(['min(s) = ',num2str(min(s))]);
    display(['s   = ',num2str(s')]);
    display(['gap = ',num2str(gap)]);
end

% no almost-stresses to blox almost-flexes
if (nw == 0 && (nv > 0 || istens))
    disp('PROBLEM: you have no stresses to block your flexes.');
    disp('Try again with different parameters.');
    return;
end


%%%==============================%%%
%%%    Find PSD stress matrix    %%%
%%%==============================%%%

% smallest singular value
if(dim == 3) niso = 3*n-6; end
if(dim == 2) niso = 2*n-3; end
sig0 = s(1 + max(0,nb-niso));


%  Case 1: no almost-flexes
if(nv == 0)
    N = null(full(C));  % orthogonal basis of space perpendicular to L
    sig0rc = min(svd(full(R*N)));  % modified July 8, 2019
    lam = 2*sig0rc^2;
    lam0 = lam;     % meaningless in this case
    kappa = 1;
    mu0 = 0;
    
    if(isbar)
        w = zeros(nb,1);
    end
    
    if(istens)
        [X,t,alph,Mevals] = pss_maxeig_tensegrity(V,W,types,dim,a,ap,ifquiet);
        
        % check if it worked
        if(isnan(alph))
            disp('FAILED to find a pss stress for this tensegrity, nv=0.');
            disp('Try again with different parameters.');
            return;
        end
        % construct optimal stress
        w = W*alph';
    end
end


% Case 2: Have almost-stresses and almost-flexes
if (nw > 0 && nv > 0 )
    
    % Do optimization to find maximally pss matrix
    if(isbar)
        [X,t,alph,Mevals] = pss_maxeig(V,W,dim,a,ap,ifquiet);
    end
    if(istens)
        [X,t,alph,Mevals] = pss_maxeig_tensegrity(V,W,types,dim,a,ap,ifquiet);
    end
    
    % Check that it was successful
    toleig = 1e-12;  % tolerance for minimum eigenvalue
    if(nv > 0 && (isnan(t) || t < toleig))  % not sure what exact tolerance should be; using tol
        disp('FAILED to find a pss stress, nv>0.');
        disp('Try again with different parameters.');
        return;
    end
    
    % Construct optimal stress matrix and extract constants
    w = W*alph';
    M = stressmatrix(w,a,ap,dim);
    mu0 = min(eig(M));
    lam0 = t;
    lam = params.lamfac * lam0;
    
    % Find kappa
    kappa = pss_kappa(R,M,C,lam,ifquiet);
    
    if(ifdisp)
        disp('---  Stress Information:  ---');
        display(['alph = ',num2str(alph)]);
        if(nv>0) display(['eig(X) = ',num2str(sort(eig(X))')]); end
        display(['eig(M) = ',num2str(sort(unique(eig(M)))')]);
    end
end


%%%========================%%%
%%%   Compute constants    %%%
%%%========================%%%

mubar = 1-mu0/lam;
L = sqrt(lam / (8*z*kappa) );
pnorm = norm(x);
eta1 = 4 * norm(w'*R) / lam;
y = eta1/L * (sqrt(mubar) + eta1/L);
y2max = 1/4 + 1/8*(sqrt(mubar)*sqrt(mubar+4) - mubar);
eta2a = L/2 * (sqrt(mubar+4) - sqrt(mubar)) - eta1;
eta2b =  L/2 * (sqrt(mubar+2) - sqrt(mubar));
y3max = 1/9 * (8/3 + sqrt(mubar)*sqrt(mubar + 8/3) - mubar);
eta3 = L/2 * (sqrt(mubar + 8/3) - sqrt(mubar));
if (kappa > 0)
    femin = @(x) sqrt( (norm(w)/kappa)^2 + (2/3)*lam/kappa * x.^2.*(1-(3/2)*eta1./x) )...
        - norm(w) / kappa;
else
    femin = @(x) 1/3 * lam / norm(w) * x.^2.*(1-(3/2)*eta1./x);
end
emin3 = femin(eta3);
ybad = 1/2 * eta1/L * (pnorm/L + eta1/L);
ypss =  eta1/L * (sqrt(mubar) + (3/2)*eta1/L + (1/2)*pnorm/L);

fetamax = @(x) (3/4)*eta1 + sqrt( ((3/4)*eta1)^2 + 3*(norm(w)/lam)*x + (3/2)*kappa/lam * x.^2 );


%%%==========================%%%
%%%   Display information    %%%
%%%==========================%%%

if(ifdisp)
    disp('---  Fundamentals:  ---');
    disp(['mubar = ',num2str(mubar)]);
    disp(['L     = ',num2str(L)]);
    disp(['|p|   = ',num2str(pnorm)]);
    disp(['kappa = ',num2str(kappa)]);
    disp(['lam0  = ',num2str(lam0)]);
    disp(['mu0   = ',num2str(mu0)]);
    disp('---  Constants for theorems to hold:  ---');
    disp(['D                 = ',num2str(y)]);
    disp(['D1max (Theorem 1) = 0.5']);
    disp(['D2max (Theorem 2) = ',num2str(y2max)]);
    disp(['D3max (Theorem 3) = ',num2str(y3max)]);
    disp(['D3max (Theorem 4) = ',num2str(ypss)]);
    disp(['Dbad  (easy bound)= ',num2str(ybad)]);
    
    
    disp('---  Radii:  ---');
    disp(['eta1  = ',num2str(eta1)]);
    disp(['eta2a = ',num2str(eta2a)]);
    disp(['eta2b = ',num2str(eta2b)]);
    disp(['eta3  = ',num2str(eta3)]);
    disp('---  Edge norm:  ---');
    disp(['emin3    = ',num2str(emin3)]);
    disp(['min squared edge = ',num2str(min(2*abs(w(abs(w)>1e-7))/kappa))]);
end




%%%========================%%%
%%%       Plot             %%%
%%%========================%%%

% Plot, if desired
if(ifplot)
    
    Vpl = V  ;
    
    if(whichWplot == 1)
        Wpl = W;
    elseif(whichWplot == 2)
        Wpl = W*alph';
    else
        Wpl = [];
    end
    
    % which flex(es) to plot
    if(~isempty(Vpl) )
        if(~isscalar(ivplot) || (ivplot ~= inf && ivplot ~=0))
            vplot = Vpl(:,ivplot);
        elseif(ivplot == inf)
            vplot = Vpl;
        elseif(ivplot == 0)
            vplot = [];
        else
            vplot = [];
        end
    else
        vplot = [];
    end
    % which stress to plot
    if(~isempty(Wpl))
        if(iwplot ~= 0)
            wplot = Wpl(:,iwplot);
        else
            wplot = [];
        end
    else
        wplot = [];
    end
    
    % plot
    whichfig = 1;
    if(dim == 2)
        tic
        plot_stressflex_2d(x,a,vplot,wplot,davg,whichfig,ap,lattice);
        toc
    end
    if(dim == 3)
        opts = struct('xth',0,'yth',0,'fig',whichfig);  % add options here (see plotcluster2b)
        plot_stressflex_3d(x,a,vplot,wplot,opts);
    end
end



