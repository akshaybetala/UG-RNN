function f=power(f,alpha)
% power : compute the factor defined by an elementwise power
% f = power(F,a) computes the factor defined by f(x)=F(x)^a 

% (c) Alexander Ihler 2010

if (~isscalar(alpha) || ~isnumeric(alpha)) error('Unknown exponent type'); end;

f.t = real(exp( log(f.t) .* alpha));
%f.t = f.t .^ alpha;

