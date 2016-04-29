function f=log(f)
% log(f): elementwise natural logarithm
% flog = log(forig) 
%   Compute a new function of the same variables, whos value is the natural-log of the original function
%

% (c) Alexander Ihler 2010

f.v = f.v+0;
f.t = log(f.t);

