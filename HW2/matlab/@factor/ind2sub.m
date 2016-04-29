function idx = ind2sub(F,ind)
% ind2sub: convert a linear (scalar) index to a subscript (cell array of indices)
% sub = ind2sub(F,ind) : compute a subscript (or variable configuration) for a given linear index 'ind'
%   As with the built-in, returns the answer in the form of a cell array
%   Uses matlab ind2sub so may not be the fastest
% See also: ind2subv, subv2ind

% (c) Alexander Ihler 2010

% TODO: vector/matrix of indices
idx=cell(1,length(F.v)); 
[idx{:}]=ind2sub(size(F.t),ind); 
%sub = [idx{:}];

