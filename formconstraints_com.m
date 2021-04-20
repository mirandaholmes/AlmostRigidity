% form center of mass constraints

% created jan 25 2018

function L = formconstraints_com(x,dim)

n = length(x)/dim;
L = sparse(repmat([1:dim]',n,1),1:n*dim,ones(size(x)),dim,n*dim);  

