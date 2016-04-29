function f=zeros(f)
% create a factor of all zeros
% f = zeros(F) : create a new factor over the same variables as F, but with f(x)=0 for all x
 f.v = f.v+0;
 f.t = zeros(size(f.t));

