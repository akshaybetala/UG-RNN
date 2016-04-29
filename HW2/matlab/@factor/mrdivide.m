function f = mrdivide(f1,f2)
% mrdivide: compute factor which is the ratio of two factors (zeros get special treatment)
% FIX : fnew = rdivide(f,g) : return a factor which is the ratio of two factors f,g
% !!! mrdivide: 0/x = 0 even if x is +/- infinity, nan, etc.
% See also: plus, minus, times

% (c) Alexander Ihler 2010


  if (~isa(f1,'factor')||~isa(f2,'factor'))                  % check for scalar products
    if (isnumeric(f1)&&isscalar(f1)) f=f2; t=f.t; t=f1./t; t(isnan(t))=0; f.t=t; return; end;
    if (isnumeric(f2)&&isscalar(f2)) f=f1; t=f.t; t=t./f2; t(isnan(t))=0; f.t=t; return; end;
    error('Unknown input type!');                                           % something else; error
  end;
  v1=f1.v; v2=f2.v; t1=f1.t; t2=f2.t;
  if (isscalar(t1)) f=f2; f.t=t1 ./ t2; f.t(isnan(f.t))=0; return; end;     % or a scalar factor class
  if (isscalar(t2)) f=f1; f.t=t1 ./ t2; f.t(isnan(f.t))=0; return; end;     % 

  f=f1; 
  [v which1 which2] = vunionmemb(v1,v2);                 % otherwise, get the new variables & sizes:
  size1=size(t1); size2=size(t2);                       % and get their variable's sizes
  if (isscalar(v1)) size1=size1(1); end;                % if there's only one variable, be careful
  if (isscalar(v2)) size2=size2(1); end;
  sizes=[which1 1]; sizes(which1)=size1; sizes(which2)=size2;         % wherever they are in the new list, set their sizes
  s11=sizes; s21=sizes; s11(~which1)=1; s21(~which2)=1;
  t=bsxfun(@rdivide,reshape(t1,s11),reshape(t2,s21));
  t(isnan(t))=0;
  f.t=t; f.v=v;


