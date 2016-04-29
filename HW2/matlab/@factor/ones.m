function f=ones(f)
% return a new factor of all ones
% f=ones(F) : create a new factor over the same variables as F, but with f(x)=1 for all x
 f.v = f.v+0;
 f.t = ones(size(f.t));
