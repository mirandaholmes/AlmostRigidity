% Make a large matrix with given matrix repeated on diagonal

% Created Jan 27 2011


function Rd = bdiag(R,n)

dum=cell(1,n);
[dum{:}]=deal(R);
Rd = blkdiag(dum{:});