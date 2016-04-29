function f = mtimes(f1,f2)
% mtimes: compute factor which is the product of two factors (zeros get special treatment)
% FIX : fnew = times(f,g) : return a factor which is the product of two factors f,g
% !!! mtimes: defines 0*x = 0  even if x is +/- infinity, NaN, etc.
% See also: plus, minus, rdivide

% (c) Alexander Ihler 2010

  if (~isa(f1,'factor')||~isa(f2,'factor'))                  % check for scalar products
    if (isnumeric(f1) && isscalar(f1)) f=f2; t=f.t; t=t.*f1; t(isnan(t))=0; f.t=t; return; end;  
    if (isnumeric(f2) && isscalar(f2)) f=f1; t=f.t; t=t.*f2; t(isnan(t))=0; f.t=t; return; end;  %
  error('Unknown input type');                       % else something wierd: error
  end;
  v1=f1.v; v2=f2.v; t1=f1.t; t2=f2.t;
  if (isscalar(t1)) f=f2; t1=t1.*t2; t1(isnan(t1))=0; f.t=t1; return; end;  % or a scalar factor class
  if (isscalar(t2)) f=f1; t1=t1.*t2; t1(isnan(t1))=0; f.t=t1; return; end;  % 

  f=f1;
  [v which1 which2] = vunionmemb(v1,v2);        % otherwise, get the full set of vars, their membership,
%  v=vunion(v1,v2);                 % otherwise, get the new variables & sizes:
%  which1 = vmember(v,v1);              % figure out which variable ids came from which
%  which2 = vmember(v,v2);              %   input factor
  size1=size(t1); size2=size(t2);          % and get their variable's sizes
  if (isscalar(v1)) size1=size1(1); end;        % if there's only one variable, be careful
  if (isscalar(v2)) size2=size2(1); end;
  sizes=[which1 1]; sizes(which1)=size1; sizes(which2)=size2;  % wherever they are in the new lis, set their sizes
  s11=sizes; s21=sizes; s11(~which1)=1; s21(~which2)=1;
  t=bsxfun(@times,reshape(t1,s11),reshape(t2,s21));
  t(isnan(t))=0;
  f.t=t; f.v=v;


function c=mytimes(a,b)
  c=a.*b;
  c(a==0|b==0)=0;
