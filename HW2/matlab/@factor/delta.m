function f=delta(F,Xv)
% delta : create appropriately sized Kronecker delta function factor
% f=delta(F,tuple) returns a factor on the same variables as F, but with
%   f(tuple)=1 and f(x)=0 otherwise.
%

t=zeros(size(F.t));
Xi=num2cell( Xv );
%Xi=num2cell( ind2sub(size(F.t),Xv) );
t(Xi{:})=1;
f=F;
f.v=f.v+0;
f.t=t;
%f=factor(F.v+0,t);

