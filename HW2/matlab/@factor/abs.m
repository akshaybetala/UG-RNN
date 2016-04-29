function f=abs(f)
% abs(F): elementwise absolute value
% f = abs(F) 
%   Compute a new function of the same variables, whos value is the absolute value of the original function
%

% (c) Alexander Ihler 2011

f.v = f.v+0;
f.t = abs(f.t);

