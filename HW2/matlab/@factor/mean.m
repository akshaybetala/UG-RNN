function f=mean(varargin)
% mean: compute the arithmetic mean of a collection of factors
% f = mean(f1,...,fn) computes the mean function f(x) = (1/n) \sum_{i=1}^n fi(x)
% for a collection of factors f1...fn
%

% (c) Alexander Ihler 2010

 factors=varargin;
 Nf=length(factors);
 f=factors{1};
 for i=2:Nf, f=f+factors{i}; end;
 f=f./Nf;

