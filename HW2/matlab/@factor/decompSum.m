function Flist = decompose(F, vlist, method)
% approximately decompose a factor into a product of smaller factors
% FList = decompSum(F, varList, method) : use one of several methods to approximate F as \sum_i f_i(v_i)
%   varList: a cell array of variable sets (uint32) over which each f_i are defined
%   FList  : a cell array of output factors
%   method can be one of: 
%      'L2'     : min L2 norm on log(F); 
%      'L2+hpm' : "" + scalar correction for Linf norm
%      'L2+mas' : "" + scalar (exponent) correction for MAS error

nF = length(vlist);

switch(lower(method))

case 'l2orig',
  v = variables(F); d = dims(F);
  vorth = vlist; vdone = uint32([]);
  Flist = cell(1,nF);
  C = (nF-1)*sum(F,v)/nF/prod(d);
  for j=1:nF, vorth{j}=vdiff(vlist{j},vdone); vdone=vunion(vdone,vlist{j}); end;
  for j=1:nF, 
    mem = vmember(v,vorth{j}); D=prod(d(~mem));
    Flist{j} = marginal(F,vorth{j}) / D - C;
  end;

case 'l2',
  v=variables(F); d=dims(F); Flist=cell(1,nF);
  Cn = sum(F,v); Cd=prod(d);
  %C = (nF-1)*sum(logF,v)/nF/prod(d);
  for j=1:nF,
    mem = vmember(v,vlist{j}); D=prod(d(~mem));
    Flist{j} = marginal(F,vlist{j}) / D - Cn/Cd*(1-1/(nF-j+1));
    %tmp2 = marginal(F,vdiff(v,vlist{j})) / prod(d(mem)) - C;
    F = F - Flist{j};
    Cn=Cn-sum(Flist{j},vlist{j})*D;
  end;

case 'l2+hpm',
    Flist = decompSum(F,vlist,'l2');
    F2=Flist{1}; for j=2:nF, F2=F2+Flist{j}; end;
    tmp = table( F-F2 ); mx=max(tmp(:)); mn=min(tmp(:));
    for j=1:nF, Flist{j}=Flist{j} + (mx+mn)/2/nF; end;

case 'l2+mas',
    Flist = decompSum(F,vlist,'l2');
    F2=Flist{1}; for j=2:nF, F2=F2+Flist{j}; end;
    % correct with an exponent
    tmp = table(log( F2/F )); mx=max(tmp(:)); mn=min(tmp(:));
    for j=1:nF, Flist{j}=Flist{j} * exp(-(mx+mn)/2/nF); end;


otherwise, error('Decomposition method not implemented.');
end;


