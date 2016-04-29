function f=geomean(varargin)
% Geometric mean of a collection of factors
% f=geomean(f1,...,fn) : compute the geometric mean of a collection of factors f1..fn
%

% (c) Alexander Ihler 2010

factors=varargin;
 Nf=length(factors);
 f=log(factors{1}); 
 for i=2:Nf, f=f+log(factors{i}); end;
 f=f./Nf; 
 f=exp(f); 

