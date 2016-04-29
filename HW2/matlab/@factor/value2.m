function v=value2(F,vars,tuples)
% return the value of F(x) evaluated at a set of configurations x
% v=value2(F,x,tuples) : return F(tuple) for each vector of tuples 
%    tuples can be Nt x Nv, where Nv=nvar(F)
%    (if tuples is Nt x 1, converts from linear index to subscript form first; see ind2sub)
% TODO: comment, compare to value

% (c) Alexander Ihler 2010, 2014

[nil, inF, inV] = vunionmemb(variables(F),vars);
%vars, variables(F), inF, inV, 
if (~all(inV)) error('vars should contain all variables in F!'); end;
tuples = tuples(:,inF);
% TODO: relax sorted requirement?


Nt=size(tuples,1); v=zeros(Nt,1); 
Nv=size(tuples,2);
assert(Nv==nvar(F) || Nv==1,'Value: configuration''s number of columns must equal one or the number of variables');

if (Nt < 200) 
  tups = num2cell(tuples);
  for i=1:Nt
  %  tup=num2cell(tuples(i,:));
  %  v(i)=F.t(tup{:});
    v(i)=F.t(tups{i,:});
  end;
else 
  for v=1:Nv,
    tups{v}=tuples(:,v);
  end;
  idx = sub2ind(size(F.t),tups{:});
  v = F.t(idx);
end

