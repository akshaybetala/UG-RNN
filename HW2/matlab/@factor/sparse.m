function sf = sparse(f)
% sparse : Convert factor to sparse spfactor class

% (c) Alexander Ihler 2011

tuples = cell(1,numel(f)); values=f.t(:)'; 
for i=1:numel(f),
  tuples{i}=ind2subv(f,i);
end;
sf = spfactor(f.v,size(f.t),0,tuples,values);

