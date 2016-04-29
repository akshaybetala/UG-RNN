function s = sample(f,nSamples)
% sample : draw a random sample from a factor, treated as a probability distribution
% s = sample(f,nSamples) : draw nSamples samples from the distribution defined by the (normalized) factor f

% TODO: slow
% TODO: add conditional sampling support?
% TODO: all-zeros factor sampling?

% (c) Alexander Ihler 2005, 2012

if (nargin < 2) nSamples = 1; end;
%s = zeros(nSamples, length(f.v));  % (???) s = zeros(length(f.v),nSamples);
%dim = dims(f);
cdf = [0;cumsum(f.t(:))]; 
[t2,idx] = histc(rand(1,nSamples)*cdf(end),[cdf(1:end-1);inf]);
s = ind2subv(f,idx);
%subs=cell(1,length(dim));
%[subs{:}]=ind2sub( dim , idx );
%s = reshape([subs{:}],nSamples,length(dim));

% 
% % (1) Should sort for large enough nSamples
% % (2) Use ind2sub or write something faster
% tmp=tmp/tmp(end);
% r = rand(1,nSamples);
% mul0 = prod(dim(1:end-1));
% for i=1:nSamples,
  % ind = find(tmp>=r(i),1,'first')-1;
  % mul = mul0;
  % for j=length(dim):-1:1
    % s(j,i) = fix(ind/mul)+1; 
    % ind = ind - (s(j,i)-1)*mul;
    % if (j~=1) mul = mul / dim(j-1); end;
  % end;    
% end;
% s=s';

