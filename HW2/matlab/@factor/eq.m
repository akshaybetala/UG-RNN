function b = eq(f1,f2,tol)
%  f==g        : test for equality of two factors (or a factor and scalar)
%  eq(f,g,tol) : test for equality up to given tolerance
  if (nargin == 3)
    t=table(abs(f1-f2))>=tol; b=any(t(:));
    return;
  end;

  if (~isa(f1,'factor')||~isa(f2,'factor'))                                 % check for scalar products
    if (isnumeric(f1) && isscalar(f1)) t=f2.t; s=f1; b=~(any(t~=s)); return; end;
    if (isnumeric(f2) && isscalar(f2)) t=f1.t; s=f2; b=~(any(t~=s)); return; end;
    error('Unknown input type');                                            % else something wierd: error
  else 
    t=f1.t; s=f2.t;
    if (any(size(t)~=size(s))) b=false;
    elseif (any(t~=s)) b=false;
    else b=true;
    end;
  end;

