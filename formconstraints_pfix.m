% form linear constraints corresponding to fixing particles

% created jan 25 2018

function L = formconstraints_pfix(x,pfix,dim)

if(dim == 2); nc = 3; end
if(dim == 3); nc = 6; end
    
L = sparse(nc,length(x));  % constraint matrix

if(dim == 3)
    p1 = pfix(1);
    p2 = pfix(2);
    p3 = pfix(3);
    L(1:3,3*p1-2:3*p1) = eye(3);  % fix xyz of p1
    L(4:5,3*p2-1:3*p2) = eye(2);  % fix yz of p2
    L(6,3*p3) = 1;                % fix z of p3
end

if(dim == 2)
    p1 = pfix(1);
    p2 = pfix(2);
    L(1:2,2*p1-1:2*p1) = eye(2);  % fix xy of p1
    L(3,2*p2) = 1;                % fix y of p2
end

