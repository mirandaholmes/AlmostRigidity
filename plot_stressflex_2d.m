% plot framework and given flex, stress, in 2d

% created jan 26 2018


function plot_stressflex_2d(points,a,vpl,wpl,sv,whichfig,varargin)

n = size(a,1);

if(nargin > 6)
    ap = varargin{1};
    lattice = varargin{2};
    if(~isempty(lattice))
        u1 = lattice(:,1);
        u2 = lattice(:,2);
    end
else
    ap = [];
    lattice = [];
end

% color map for stress 0.2 for final plot, 0.7 for now
lcolr = 0.7*[1 1 1];    % default colour if stress is below threshhold
lstyl = '-';
lw = 3;
    cpos  =     [0 0 1];       % color for positive stress 
    cneg  =     [1 0 0];       % color for negative stress
czero = 1*[1 1 1];     % colour for zero stress


wtol = 0.1;%0.02;  % threshhold at which to plot stress (0.02 for jayson data)
%sv = 1;     % size of motion for v
msize = 8;   % marker size for points
iftext = 1;
if(n > 200)
    lw = 2;
    msize = 6;
end
if(n > 400)
    msize = 2;
    iftext = 0;
    lw = 1;
end


cmap = lines(7);  % default color map
cblue1 = [0    0.4470    0.7410];
cblue2 = [0.6    0.9    1];
cblue3 = [0.2    0.9    0.7];
corange = [0.8500    0.3250    0.0980];
cred = [0.6350    0.0780    0.1840];
cpurple = [0.5,0,0.5];

 % color map for velocity
cmapv = {'green','yellow',corange,'blue',cpurple,'red'}; 
ncmap = length(cmapv);
nv = size(vpl,2);
nw = size(wpl,2);




% -----   Construct needed variables   -----

x = points(1:2:end);
y = points(2:2:end);
[rr,cc] = find(triu(a));
[rp,cp,sp] = find(ap);
edges = [rr,cc];
edges_p = [rp,cp];
edges_all = [edges;edges_p];
nb = length(rr);  % number of regular edges
nbp = length(rp);    % number of periodic edges


% -----  Set up plot  -----

figure(whichfig);
clf
hold on

% -----   Plot edges and stress  -----
if(nw == 0)
    wline = zeros(nb+nbp,1);  
else
    wline = wpl / max(abs(wpl));  % scale so maximum intensity is 1
end



% Regular edges
i1 = edges_all(1:nb,1);
i2 = edges_all(1:nb,2);
x1 = x(i1); x2 = x(i2);
y1 = y(i1); y2 = y(i2);
ws = wline(1:nb);  % value of stress on this bond
lcolrw = repmat(lcolr,nb,1);
for ie=1:nb
    if(ws(ie) > wtol)
        lcolrw(ie,:) = ws(ie)*cpos + (1-ws(ie))*czero;
    elseif(ws(ie) < -wtol)
        lcolrw(ie,:) = (-ws(ie))*cneg + (1+ws(ie))*czero;
    end
end
h = plot([x1,x2]',[y1,y2]','Linewidth',lw,'LineStyle',lstyl);
set(h,{'color'},num2cell(lcolrw,2));



% Periodic edges
copies_x = [];  % keeps track of periodic copies, for plotting  
copies_y = [];
for ie=(nb+1):nb+nbp
    i1 = edges_all(ie,1);
    i2 = edges_all(ie,2);
    x1 = x(i1); x2 = x(i2);
    y1 = y(i1); y2 = y(i2);
    
    % get colour of edge, based on value of stress
    ws = wline(ie);  % value of stress on this bond
    if(abs(ws) > wtol)
        if(ws >= 0)
            lcolrw = ws*cpos + (1-ws)*czero;
        else
            lcolrw = (-ws)*cneg + (1+ws)*czero;
        end
    else
        lcolrw = lcolr;
    end
    
    % shift endpoints, if periodic
    if(ie > nb)
        if(sp(ie-nb) < 0)
            x1t = x1 + u1(1);
            y1t = y1 + u1(2);
            x2t = x2 - u1(1);
            y2t = y2 - u1(2);
        else
            x1t = x1 + u2(1);
            y1t = y1 + u2(2);
            x2t = x2 - u2(1);
            y2t = y2 - u2(2); 
        end
        copies_x = [copies_x;x1t;x2t];
        copies_y = [copies_y;y1t;y2t];
    end
    

    % draw the lines
    if(ie <= nb)  % regular bond
       line([x1,x2],[y1,y2],'color',lcolrw,'Linewidth',lw,'LineStyle',lstyl);
    else  % periodic bond: draw copies coming from both ends
       line([x1t,x2],[y1t,y2],'color',lcolrw,'Linewidth',lw,'LineStyle',lstyl);
       line([x1,x2t],[y1,y2t],'color',lcolrw,'Linewidth',lw,'LineStyle',lstyl);
    end
end


% -----   Plot flex -----
if(nv > 0)
    scale = 0;  % don't automatically rescale the vectors
    for jv=1:nv
        vv = vpl(:,jv);
        vv2 = reshape(vv,2,n);
        normv = sum(vv2.^2,1);
        vv = vv/sqrt(max(normv))*sv;  % make largest vector have given length
        quiver(x,y,vv(1:2:end),vv(2:2:end),scale,...
            'color',cmapv{mod(jv-1,ncmap)+1},'MaxHeadSize',2,'LineWidth',1.5);
    end
end



% -----   Plot points   -----
% Plot interior points
plot(x,y,'o','Markersize',msize,'MarkerEdgeColor','k','MarkerFaceColor',cblue1);

% Plot periodic copies
plot(copies_x(1:2:end),copies_y(1:2:end),'o','Markersize',msize,'MarkerEdgeColor','k','MarkerFaceColor',cblue2);
plot(copies_x(2:2:end),copies_y(2:2:end),'o','Markersize',msize,'MarkerEdgeColor','k','MarkerFaceColor',cblue3);





% -----   Plot text   -----
if(iftext)
    tshift = 0.04;
    for jx=1:length(x)
        text(x(jx)+tshift,y(jx),num2str(jx),'HorizontalAlignment','left')
    end
end

hold off
axis equal


