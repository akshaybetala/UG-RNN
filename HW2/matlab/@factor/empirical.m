function F = empirical(F,data)
% empirical(F, data) : construct an empirical distribution of the data for the variables in F
%   data should be (# samples) x (# variables) and take values 1..Di where D=dims(F)
% 
% Example:  Fhat = empirical(F, sample(F,1000)); 

% (c) Alexander Ihler 2012

if (isempty(data)) F=factor; return; end;
if (nvar(F) ~= size(data,2)) error('empirical: data must have the same number of columns as F has variables'); end;

F.t = zeros([dims(F) 1 1]);
if (numel(F) < size(data,1))   % if fewer configurations than samples,
  for i=1:numel(F),            %  count by configurations (?)
    Xi = ind2subv(F,i);
    F.t(i) = sum( all(bsxfun(@eq,data,Xi),2) );
  end;
else
  for i=1:size(data,1)         % else walk through the data and convert to index
    idx = subv2ind(F,data(i,:));
    F.t(idx) = F.t(idx) + 1;
  end;
end;
F.t=F.t/size(data,1);

