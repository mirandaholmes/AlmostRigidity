% plot framework and given flex, stress, in 3d

% created jan 26 2018

% TO DO: 
% -- incorporate periodic constraints

function plot_stressflex_3d(x,a,vpl,wpl,varargin)

sv = 2;  % size of arrows for flex

corange = [0.8500    0.3250    0.0980];
cpurple = 0.9*[1 0 1];
cmapv = {'green','yellow','blue',corange,cpurple};  % color map for velocity
ncmap = length(cmapv);

nb = sum(sum(triu(a)));   % # of bonds
n = length(x)/3;          % # of particles
nv = size(vpl,2);         % # of flexes to plot

% colour for spheres (grey)
scolr2 = 0.5*[1 1 1];  

% colour for bars, based on stress
if(isempty(wpl))
    lcolr = 0.7*[1,1,1];  % default; grey
else
    lcolr = zeros(nb,3);   % coloured bars
    wline = wpl / max(abs(wpl));  % scales so maximum intensity is 1
    cpos  =     [0 0 1];       % color for positive stress 
    cneg  =     [1 0 0];       % color for negative stress
    czero = 1*[1 1 1];  % color for zero stress
    for jb=1:nb
        ws = wline(jb);
        if(ws >= 0)
            lcolr(jb,:) = ws*cpos + (1-ws)*czero;
        else
            lcolr(jb,:) = (-ws)*cneg + (1+ws)*czero;
        end
    end
end


% Options for plotting
 opts1 = struct('srad',0.1,'lcolr',lcolr,'lrad',0.05,...
                'scolr',scolr2,'salph',1,'ifgrid',0,'iftext',1,...
                'a',a);
            
            
 % Extract user-defined options
if(nargin > 4)
    % Convert varargin to struct, if it's not already (previous versions it
    % used to be a struct, now it's a cell array)
    varargin = varargin{1};  % hack for my laptop
    if(~isstruct(varargin))
        opts2=cell2struct(varargin(2:2:end), varargin(1:2:end), 2);
    else
        opts2 = varargin;%varargin{1};
    end
    
    % Merge options, using user-defined as the last ones so they get priority
    opts = catstruct(opts1,opts2);
else
    opts = opts1;
end


% Plot the cluster
 [h2,pos2,x2] = plotcluster2b(x,opts);
 
 
 % Rotate arrows, if needed
 if(isfield(opts,'zth') && ~isempty(opts.zth))
     zth = opts.zth;
 else
     zth = 0;
 end
 if(isfield(opts,'yth') && ~isempty(opts.yth))
     yth = opts.yth;
 else 
     yth = 0;
 end
 if(isfield(opts,'xth') && ~isempty(opts.xth))
     xth = opts.xth;
 else 
     xth = 0;
 end
 
% Rotation matrices
Rx = @(th) ([1, 0, 0, ; ...
             0, cos(th), -sin(th) ;...
             0, sin(th), cos(th) ]);
Ry = @(th) ([cos(th), 0, sin(th);...
             0, 1, 0;...
            -sin(th), 0, cos(th)]);
Rz = @(th) ([cos(th), -sin(th), 0; ...
             sin(th), cos(th), 0; ...
             0, 0, 1]);    
         
 
% Now add arrows
for jv=1:nv
    vv = vpl(:,jv);   % jv'th direction
    V = Rz(zth)*Ry(yth)*Rx(xth)*reshape(vv,3,n);  % rotate vector
    vout = reshape(V,3*n,1);
    for jx=1:n
        p0  = x2(3*jx-2:3*jx);    % start point (need x2 = output of plotcluster2b)
        dir = vout(3*jx-2:3*jx);    % direction
        p1  = p0 + sv*dir;         % end point
        mArrow3(p0,p1,'color',cmapv{mod(jv-1,ncmap)+1});%,'stemWidth',0.015,'tipWidth',0.06);
    end
end

 
 
 
 


