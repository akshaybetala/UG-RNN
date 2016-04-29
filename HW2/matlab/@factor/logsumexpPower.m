function f = logsumexpPower(F,var,p)
% logsumexpPower : weighted or power log-sum-exp elimination operator
% f = logsumexpPower(F,vars,p) eliminates 'vars' in F by  f(x) = (1/p) LSE_{vars} (p*F(x))

% (c) Alexander Ihler 2010

  switch (p)
  case  1.0, f=logsumexp(F,var); 
  case  inf, f=max(F,var); 
  case -inf, f=min(F,var); 
  otherwise, f=logsumexp(F*p,var)/p; 
  end;

