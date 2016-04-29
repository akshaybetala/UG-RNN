function Contents(f)
% List methods availble for the factor class (help factor/methods)
% Methods available for the factor class:
% Constructors:
%   factor(vars,table)          : create a factor over variables "vars" with entries given by "table"
%   readErgo(factor,filename)   : read a list of factors from various formats.
%   readUAI10(factor,filename)
%   readWCSP(factor,filename)
%   changeVars                  : create a new factor with the same table but different arguments
%   normalize                   : create a new factor that is "normalized" in various ways
%   empirical(factor,data)      : create the empirical distribution of the data using input factor as template
%
% Accessors:
%   variables(F)                : get the list of arguments for F
%   table(F)                    : get the table of values for F
%   nvar(F)                     : number of arguments
%   dims(F)                     : dimensions of argument variables
%   value(F,tuples)             : evaluate F for the configurations "tuples" (Ntup x Nvar)
%   display(F)                  : print information about the factor
%   
% Checks and conversions:
%   ==,~=                       : test for equality of function itself
%   >,>=,<,<=                   : element-wise comparison operators
%   isempty(F)                  : false if F has values (constant or otherwise)
%   isnan(F)                    : true if any entries of F are NAN
%   isfinite(F)                 : false if any entries of F are non-finite
%   isscalar                    : true if F is a constant value (no arguments)
%   numel                       : number of elements required to specify F (=prod of dims)
%   ind2sub, ind2subv, subv2ind : convert between tuple- and index-based representations
%
% Elementwise operators
%   abs                         : take absolute value, f(x) = |F(x)|
%   exp                         : take exponential, f(x) = exp(F(x))
%   power                       : take power, f(x) = F(x)^a
%   log, log2, log10            : take logarithms, f(x) = log(F(x))
%   
% Operators
%   plus, minus, times, rdivide : standard operators on factors (+,-,.*,./)
%   mtimes, mrdivide            : (*,/) - like .*, ./ except behave differently on zero-values
%   sum, sumPower               : summation and "power" summation
%   logsumexp                   : summation for log-domain factors
%   max, min                    : max or min elimination operators
%   argmax, argmin              : find one or more optimizing values of F
%   sample                      : draw a random value from F (normalizes first)
%   condition, slice            : construct a sub-table conditioned on some variables' states
%   marginal                    : marginalize out some variables
%   maxmarginal, minmarginal    : max/minimize out some variables
%   
% Evaluators
%   distance(F1,F2)             : compute the distance between two factors under various norms
%   norm(F)                     : compute the length of F under various norms
%   entropy(F)                  : compute the entropy of F (normalizes first)
%
% Methods on collections of factors
%   mean(f1,...,fn)             :  compute F = 1/N \sum_i f_i(x)
%   geomean(f1,...,fn)          :  compute F = (prod_i f_i(x))^(1/N)
%   

  help factor/Contents







