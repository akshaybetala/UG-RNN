function F = substitute(F,a,b)
% F=substitute(F,a,b) : convert instances of value "a" to value "b" in table F

if (isnan(a))
  F.t( isnan(F.t) ) = b;
else 
  F.t( F.t==a ) = b;
end;  

