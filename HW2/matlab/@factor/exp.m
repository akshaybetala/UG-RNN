function f=exp(f)
% Element-wise exponentiation
% fexp = exp(forig) 
% Compute a new function of the same variables, whos value is the exponential of the original function
%

% (c) Alexander Ihler 2010

f.v = f.v+0;
f.t = exp(f.t);

