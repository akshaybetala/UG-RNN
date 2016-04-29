function N=norm(f,type)
% norm: compute the norm of the function/factor under various metrics
% dist = norm(f1,method)
%  Compute the norm (length) of a factor under some metric.  'method' can be one of:
%    'L1'   : L1 or integrated absolute value
%    'L2'   : L2 or integrated squared value
%    'Linf' : L-infinity or maximum absolute error
%    'hpm'  : Hilbert's projective metric norm

% (c) Alexander Ihler 2010


if (nargin==1) type='L2'; end;

switch(lower(type))
 case 'l2',  N=sum(f.t(:).^2);
 case 'l1',  N=sum(abs(f.t(:)));
 case 'linf',N=max(abs(f.t(:)));
 case 'hpm', t=log(f.t(:)); N=max(t)-min(t);
end;

