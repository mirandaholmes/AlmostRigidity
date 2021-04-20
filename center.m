% Function to put center of mass at origin, for Planar

% Created April 25, 2012

function xc = center(x,varargin)

if(nargin > 1)
    d = varargin{1};
else
    d = 3;
end

xc = x;
for jj=1:d
    xc(jj:d:end) = xc(jj:d:end) - mean(xc(jj:d:end));
end
