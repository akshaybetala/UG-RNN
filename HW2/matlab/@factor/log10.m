function f=log10(f)
% log10(f): element-wise logarithm, base 10
% flog = log10(forig) 
%   Compute a new function of the same variables, whos value is the log-base-10 of the original function
%

% (c) Alexander Ihler 2010

f.v = f.v+0;
f.t = log10(f.t);

