function sub = ind2subv(F,ind)
% ind2subv: convert a linear (scalar) index to a subscript *vector* array of indices
% sub = ind2subv(F,ind) : compute a subscript (or variable configuration) for a given linear index 'ind'
%   returns the configuration in the form of a vector (rather than a cell array)
%   Uses matlab ind2sub so may not be the fastest
% See also: ind2sub, subv2ind

% (c) Alexander Ihler 2010

% TODO: vector/matrix of indices
if (isscalar(F)) sub=1; return; end;
idx=cell(1,length(F.v)); 
[idx{:}]=ind2sub(size(F.t),ind); 
sub = reshape([idx{:}], length(ind),length(idx));

