function f=log2(f)
% log2(f): elementwise logarithm, base two
% flog = log2(forig [,zflag] ) 
%   Compute a new function of the same variables, whos value is the log-base-2 of the original function

% (c) Alexander Ihler 2010

f.v = f.v+0;
f.t = log2(f.t);

