function f = sumPower(F,var,p)
% sumPower : weighted or power sum elimination operator
% f = sumPower(F,vars,p) eliminates 'vars' in F by  f(x) = [sum_{vars} (F(x))^(p) ]^(1/p)

% (c) Alexander Ihler 2010

  switch (p)
  case 1.0,  f=sum(F,var); return;
  case inf, f=max(F,var); return;
  case -inf, f=min(F,var); return;
  otherwise
  %f=F.^(p);    % Linear domain version
  %f=sum(f,var);
  %f=f.^(1/p);
  f=log(F).*(p); % Log-domain version
  f=logsumexp(f,var);
  f=exp(f./p);
  end;

