function Flist = decompose(F, vlist, method)
% approximately decompose a factor into a product of smaller factors
% FList = decompProd(F, varList, method) : use one of several methods to approximate F as \prod_i f_i(v_i)
%   varList: a cell array of variable sets (uint32) over which each f_i are defined
%   FList  : a cell array of output factors
%   method can be one of: 
%      'L2'     : min L2 norm on log(F); 
%      'L2+hpm' : "" + scalar correction for Linf norm
%      'L2+mas' : "" + scalar (exponent) correction for MAS error

Flist = decompSum(log(F),vlist,method);
for j=1:length(Flist), Flist{j}=exp(Flist{j}); end;

