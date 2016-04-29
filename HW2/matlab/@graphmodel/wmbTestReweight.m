function gm = wmbTestReweight(gm, Task)


for i=1:length(gm.Alg.minibucket),
  nMB = length(gm.Alg.minibucket(i).nodes);
  for ni=1:nMB,
    n = gm.Alg.minibucket(i).nodes(ni);
    switch (Task),
    case 'sum+', gm.Alg.nodes(n).wt = 1.0/nMB;
    case 'sum-', if (ni==1) gm.Alg.nodes(n).wt = 2.0; else gm.Alg.nodes(n).wt = -1.0/(nMB-1); end;
    case 'max+', gm.Alg.nodes(n).wt = eps;
    end;
  end;
end;
