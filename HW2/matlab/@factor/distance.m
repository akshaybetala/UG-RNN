function N=distance(f1,f2,type)
%  Compute the distance between two factors under various distance metrics
%  N=distance(f1,f2,method), where 'method' can be one of:
%    'L1'   : L1 or integrated absolute error
%    'L2'   : L2 or integrated squared error
%    'Linf' : L-infinity or maximum absolute error
%    'KL'   : Kullback-Leibler divergence D(f1 || f2)
%    'hpm'  : Hilbert's projective metric
%

% (c) Alexander Ihler 2010

if (nargin<2) error('distance: needs at least two factor arguments'); end;
if (nargin==2) type='L2'; end;

switch(lower(type))
 case 'l2',  N=sum((f1.t(:)-f2.t(:)).^2);
 case 'l1',  N=sum(abs(f1.t(:)-f2.t(:)));
 case 'linf',N=max(abs(f1.t(:)-f2.t(:)));
 case 'kl',  N=sum(f1*log(f1/f2)); 
 %case 'kl',  t1=f1.t(:); t2=f2.t(:); t1=t1/sum(t1); t2=t2/sum(t2); N=t1'*log(t1./t2);
 %case 'kl',  t1=f1.t(:); t2=f2.t(:); t1=t1/sum(t1); t2=t2/sum(t2); N=t1'*log2(t1./t2);
 case 'hpm', t=log(f1.t(:))-log(f2.t(:)); N=max(t)-min(t);
 case 'mas', t=log(f1.t(:))./log(f2.t(:)); N=max(max(t), 1/min(t))-1;
 otherwise, error(['Unknown distance type ' type]);
end;

