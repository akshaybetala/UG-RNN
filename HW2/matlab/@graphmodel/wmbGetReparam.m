function [flist, wts] = wmbGetReparam(wmb)
% [flist wts] = wmbGetReparam(wmb)
% Extract a reparameterization of the model (list of factors) from an optimized WMB

lnZ = 0.0;
flist = {};
wts = [];

for i=1:length(wmb.Alg.order),
  for n=wmb.Alg.nodes{i},
    if (wmb.Alg.parent(n) == 0), 
      lnZ = lnZ+table(wmb.Alg.msgFwd{n});
    end;
    flist{end+1} = wmb.Alg.theta{n} - wmb.Alg.msgFwd{n};
    for c=wmb.Alg.children{n}, flist{end}=flist{end} + wmb.Alg.msgFwd{c}; end;
    wts(end+1) = wmb.Alg.wt(n);
  end;
end;

lnZ = lnZ/length(flist);
for f=1:length(flist), flist{f}=exp(flist{f}+lnZ); end;  

