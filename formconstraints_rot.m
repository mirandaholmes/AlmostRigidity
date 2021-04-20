% make constraints in the infinitesimal rotation directions

% created jan 25 2018



function L = formconstraints_rot(x,dim)

n = length(x)/dim;
th = pi/2;  

% 2d rotation matrix
R2d = ([cos(th), -sin(th); ...
        sin(th), cos(th)]);   
 
% 3d rotation matrices
b = 0;  % b = 1 for actual rotation matrices; b=0 gives derivative
Rx = ([b, 0, 0, ; ...
       0, cos(th), -sin(th) ;...
       0, sin(th), cos(th) ]);
Ry = ([cos(th), 0, sin(th);...
       0, b, 0;...
       -sin(th), 0, cos(th)]);
Rz = ([cos(th), -sin(th), 0; ...
       sin(th), cos(th), 0; ...
       0, 0, b]);  
         
% center x
if(dim == 2)
    xc = [sum(x(1:2:end));sum(x(2:2:end))]/n;
end
if(dim == 3)
    xc = [sum(x(1:3:end));sum(x(2:3:end));sum(x(3:3:end))]/n;
end
x2 = x - repmat(xc,n,1);

% infinitesimal rotations         
if(dim == 2)
    Rn = bdiag(R2d,n);
    L = (Rn*x2)';
end
if(dim == 3)
    Rxn = bdiag(Rx,n);  % Rotation (will be (Rxd*F)', etc.)
    Ryn = bdiag(Ry,n);  
    Rzn = bdiag(Rz,n);
    L = [Rxn*x2,Ryn*x2,Rzn*x2]';
end
L = sparse(L);


    
    
    