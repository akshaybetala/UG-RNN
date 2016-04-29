function ind = subv2ind(F,subv)
% convert a vector of subscript / variable values into a single (linear) index
% ind = subv2ind(F,sub) : compute a linear index for a given vector-valued variable configuration
%   Uses matlab sub2ind so may not be the fastest
% See also: ind2sub, ind2subv

% (c) Alexander Ihler 2010

ind = zeros(size(subv,1),1);
dim = size(F.t);
for i=1:size(subv,1),
 idx = num2cell(subv(i,:));
 ind(i) = sub2ind(dim,idx{:});
end;

