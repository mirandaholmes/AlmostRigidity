% Plot a cluster of spheres
% Takes many optional arguments, to modify plots
% Based on plotcluster2.m
%
% This version doesn't call any outside functions
%
% Created Nov 29, 2016
%

function [h1,pos,varargout] = plotcluster2b(x,varargin)

if(size(x,2) > 1)  % make x a vector
    x = reshape(x,size(x,1)*size(x,2));
end
n = length(x)/3;
xout = NaN(size(x));


% Colours
vlightblue = [0.0,0.5,1];
vred = [0.7,0.0,0.2];
vblue = 0.8*[0,0,1];
vgreen = [0,1,0.2];

% -----------  Set default values for all parameters  -----------

bondeps = 1e-6;
ax = 1.4*[-1 1 -1 1 -1 1];
iftight = 1;
ifcenter = 1;
iftext = 0;
iflines = 1;
ifgrid = 0;
scolr = 0.5*[1,1,1];  % color of spheres
srad = 0.25;  % radius of spheres
salph = 1;  % transparency for spheres, 0 for transparent (0.65 orig)
lalph = 1;  % transparency for lines
lightcolr = 0.5*[1 1 1];  % colour of light
lightpos = [-1 0.25 1]; %[0 0.25 1];  % light position
ambstrength = 0.5;  % intensity of ambient component of light reflected from object
specstrength = 0.8;  % intensity of specular component of reflected light
diffstrength = 1;  % intensity of diffuse component of reflected light
specexp = 2;    % specular exponent (large = small light spots)
lcolr = 0.5*[1,1,1];  % line colour (cylinder)
lrad = 0.04;  % radius of lines (cylinder)  % increased from 0.03 feb 27 2014
zth = 0;  % angle to rotate around z-axis
yth = 0;     % angle to rotate around y-axis
xth = 0;   % angle to rotate around x-axis
az = NaN; el = NaN;  % view parameters (default = 3d view)
pos = [24,10,15,15];
fig = 1;
tfs = 24;  % text font size
tshift = 0.12;  % text shift
% create text labels (default is numerical)
for ii=1:n
    textlabels{ii} = num2str(ii);
end

scolr0=[1,0,0;...
    0,1,0;...
    0,0,1;...
    1,1,0;...
    1,0,1;...
    0,1,1;...
    0.25,0.25,0;...
    0,0.25,0.25;...
    0.25,0,0.25;...
    0,0,0.5;...
    0,0.5,0;...
    0.5,0,0;...
    ];  % default color scheme

% -----------  Change defaults, if necessary  -----------

if(nargin > 1)
    % Convert varargin to struct, if it's not already (previous versions it
    % used to be a struct, now it's a cell array)
    varargin = varargin{1};  % hack for my laptop
    if(~isstruct(varargin))
        opts=cell2struct(varargin(2:2:end), varargin(1:2:end), 2);
    else
        opts = varargin;%varargin{1};
    end
    % Extract field values
    if(isfield(opts,'bondeps') && ~isempty(opts.bondeps))
        bondeps = opts.bondeps;
    end
    if(isfield(opts,'srad') && ~isempty(opts.srad))
        srad = opts.srad;
    end
    if(isfield(opts,'ax') && ~isempty(opts.ax))
        ax = opts.ax;
    end
    if(isfield(opts,'iftight') && ~isempty(opts.iftight))
        iftight = opts.iftight;
    end
    if(isfield(opts,'salph') && ~isempty(opts.salph))
        salph = opts.salph;
    end
    if(isfield(opts,'lalph') && ~isempty(opts.lalph))
        lalph = opts.lalph;
    end
    if(isfield(opts,'iftext') && ~isempty(opts.iftext))
        iftext = opts.iftext;
    end
    if(isfield(opts,'iflines') && ~isempty(opts.iflines))
        iflines = opts.iflines;
    end
    if(isfield(opts,'ifgrid') && ~isempty(opts.ifgrid))
        ifgrid  = opts.ifgrid;
    end
    if(isfield(opts,'scolr') && ~isempty(opts.scolr))
        if(opts.scolr==1)
            scolr = [scolr0;scolr0];
        else
            scolr = opts.scolr;
        end
    end
    if(isfield(opts,'lightcolr') && ~isempty(opts.lightcolr))
        lightcolr = opts.lightcolr;
    end
    if(isfield(opts,'lightpos') && ~isempty(opts.lightpos))
        lightpos = opts.lightpos;
    end
    if(isfield(opts,'lcolr') && ~isempty(opts.lcolr))
        lcolr = opts.lcolr;
    end
    if(isfield(opts,'lrad') && ~isempty(opts.lrad))
        lrad = opts.lrad;
    end
    if(isfield(opts,'zth') && ~isempty(opts.zth))
        zth = opts.zth;
    end
    if(isfield(opts,'yth') && ~isempty(opts.yth))
        yth = opts.yth;
    end
    if(isfield(opts,'xth') && ~isempty(opts.xth))
        xth = opts.xth;
    end
    if(isfield(opts,'az') && ~isempty(opts.az))
        az = opts.az;
    end
    if(isfield(opts,'el') && ~isempty(opts.el))
        el = opts.el;
    end
    if(isfield(opts,'pos') && ~isempty(opts.pos))
        pos = opts.pos;
    end
    if(isfield(opts,'fig') && ~isempty(opts.fig))
       fig = opts.fig;
    end
    if(isfield(opts,'a') && ~isempty(opts.a))
       a = opts.a;
       if(size(a,2)==1)
           a = reshape(a,n,n);
       end
    end
    if(isfield(opts,'tfs') && ~isempty(opts.tfs))
       tfs = opts.tfs;
    end
    if(isfield(opts,'tshift') && ~isempty(opts.tshift))
       tshift = opts.tshift;
    end
    if(isfield(opts,'ambstrength') && ~isempty(opts.ambstrength))
       ambstrength = opts.ambstrength;
    end
    if(isfield(opts,'specstrength') && ~isempty(opts.specstrength))
       specstrength = opts.specstrength;
    end
    if(isfield(opts,'diffstrength') && ~isempty(opts.diffstrength))
       diffstrength = opts.diffstrength;
    end
    if(isfield(opts,'specexp') && ~isempty(opts.specexp))
       specexp = opts.specexp;
    end
    if(isfield(opts,'textlabels') && ~isempty(opts.textlabels))
       textlabels = opts.textlabels;
    end
    if(isfield(opts,'ifcenter') && ~isempty(opts.ifcenter))
        ifcenter = opts.ifcenter;
    end
end

% Get adjacency matrix
if(exist('a'))
    a0 = a;
else
    %a0 = get_a_of_x(x, 1 + bondeps); 
    r = zeros(n);
    for ii=1:n
        for jj=ii+1:n
            r(ii,jj) = norm(x(3*ii-2:3*ii)-x(3*jj-2:3*jj));
            r(jj,ii) = r(ii,jj);
        end
    end
    a0 = +( (r + 2*(1+bondeps)*eye(n))< 1 + bondeps);
end


% -----------  Plot stuff  -----------

h1 = figure(fig);
clf
set(h1,'Units','centimeters');
set(h1,'Position',pos);

% Set up spheres
[sx,sy,sz]= sphere(40);  %-- sphere facets
%[Rx,Ry,Rz] = getRotations;  % rotation matrix functions
Rx = @(th) ([1, 0, 0, ; ...
             0, cos(th), -sin(th) ;...
             0, sin(th), cos(th) ]);
Ry = @(th) ([cos(th), 0, sin(th);...
             0, 1, 0;...
            -sin(th), 0, cos(th)]);
Rz = @(th) ([cos(th), -sin(th), 0; ...
             sin(th), cos(th), 0; ...
             0, 0, 1]);   
%x = center(x);
if(ifcenter == 1)
    xc = x;
    for jj=1:3
        xc(jj:3:end) = xc(jj:3:end) - mean(xc(jj:3:end));
    end
    x = xc;
end

X = Rz(zth)*Ry(yth)*Rx(xth)*reshape(x,3,n);

xout = reshape(X,3*n,1);

del = 0;
axC = axes('Position',[0+del 0+del 1-2*del 1-2*del]);  % figure goes from 0 to 1, takes full window
hold on
% plot spheres
for jn=1:n
    if(size(scolr,1)==1)
        scolt = scolr;
    else
        scolt = scolr(jn,:);
    end
    surf(sx*srad+X(1,jn), sy*srad+X(2,jn), sz*srad+X(3,jn),...
        'LineStyle','none',...
        'FaceColor',scolt,...
        'FaceAlpha',salph,...
        'DiffuseStrength',diffstrength,...
        'AmbientStrength',ambstrength,...
        'SpecularStrength',specstrength,...
        'SpecularExponent',specexp);
end

% plot lines
if(iflines)
    ind = find(triu(a0));
    siz = size(a0);
    for jb=1:length(ind)  % plot non-broken bonds
        [ir,ic] = ind2sub(siz,ind(jb));
        v1 = [X(1,ir),X(2,ir),X(3,ir)];
        v2 = [X(1,ic),X(2,ic),X(3,ic)];
        [xc,yc,zc] = cylinder2P(lrad,20,v1,v2);
        if(size(lcolr,1)==1)
            lcolt = lcolr;
        else
            lcolt = lcolr(jb,:);
        end
        surf(xc,yc,zc,...
            'LineStyle','none',...
            'FaceColor',lcolt,...
            'FaceAlpha',lalph,...
            'DiffuseStrength',diffstrength,...
            'AmbientStrength',ambstrength,...
            'SpecularStrength',specstrength,...
            'SpecularExponent',specexp);
        % Previous code
        %     line([X(1,ir) X(1,ic)],[X(2,ir) X(2,ic)],[X(3,ir) X(3,ic)],...
        %         'Color',lcolr,'LineWidth',lw);
    end
end
hold off
daspect([1,1,1]);
axis(ax);
view(3);
if(~isnan(az))
    view([az,el]);
end
light('Position',lightpos,'Style','infinit','Color',lightcolr);
lighting phong
axis off
if(ifgrid)
    axis on
    grid on
    xlabel('x');ylabel('y');zlabel('z');
else
    axes(axC);
end
if(iftight)
    axis tight
end

% sphere labels
if(iftext)
    for jn=1:n
        yourtext = textlabels(jn);  
        text(X(1,jn)-tshift,X(2,jn)-tshift,X(3,jn)+tshift,yourtext,'FontSize',tfs,...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
end

if(nargout > 2)
    varargout{1} = xout;
end

%set(h1,'PaperPosition',pos);    
    





%  CYLINDER:  A function to draw a N-sided cylinder based on the
%             generator curve in the vector R.
%
%  Usage:      [X, Y, Z] = cylinder(R, N)
%
%  Arguments:  R - The vector of radii used to define the radius of
%                  the different segments of the cylinder.
%              N - The number of points around the circumference.
%
%  Returns:    X - The x-coordinates of each facet in the cylinder.
%              Y - The y-coordinates of each facet in the cylinder.
%              Z - The z-coordinates of each facet in the cylinder.
%
%  Author:     Luigi Barone
%  Date:       9 September 2001
%  Modified:   Per Sundqvist July 2004

function [X, Y, Z] = cylinder2P(R, N,r1,r2)

% The parametric surface will consist of a series of N-sided
% polygons with successive radii given by the array R.
% Z increases in equal sized steps from 0 to 1.

% Set up an array of angles for the polygon.
theta = linspace(0,2*pi,N);

m = length(R);                 % Number of radius values
% supplied.

if m == 1                      % Only one radius value supplied.
    R = [R; R];                % Add a duplicate radius to make
    m = 2;                     % a cylinder.
end


X = zeros(m, N);             % Preallocate memory.
Y = zeros(m, N);
Z = zeros(m, N);

v=(r2-r1)/sqrt((r2-r1)*(r2-r1)');    %Normalized vector;
%cylinder axis described by: r(t)=r1+v*t for 0<t<1
R2=rand(1,3);              %linear independent vector (of v)
x2=v-R2/(R2*v');    %orthogonal vector to v
x2=x2/sqrt(x2*x2');     %orthonormal vector to v
x3=cross(v,x2);     %vector orthonormal to v and x2
x3=x3/sqrt(x3*x3');

r1x=r1(1);r1y=r1(2);r1z=r1(3);
r2x=r2(1);r2y=r2(2);r2z=r2(3);
vx=v(1);vy=v(2);vz=v(3);
x2x=x2(1);x2y=x2(2);x2z=x2(3);
x3x=x3(1);x3y=x3(2);x3z=x3(3);

time=linspace(0,1,m);
for j = 1 : m
    t=time(j);
    X(j, :) = r1x+(r2x-r1x)*t+R(j)*cos(theta)*x2x+R(j)*sin(theta)*x3x;
    Y(j, :) = r1y+(r2y-r1y)*t+R(j)*cos(theta)*x2y+R(j)*sin(theta)*x3y;
    Z(j, :) = r1z+(r2z-r1z)*t+R(j)*cos(theta)*x2z+R(j)*sin(theta)*x3z;
end

%surf(X, Y, Z);



