function gm = wmbSetBucketWeight(gm, bucket, weight)
% set (mini)bucket weights.  weight can be a vector (# of minibuckets) or scalar (+1,eps,-1)=>evenly distributed 
% TODO: improve comments

  nMB = length(gm.Alg.minibucket(bucket).nodes);
  if (length(weight)~=1 && length(weight)~=nMB) error('Wrong size weight vector?'); end;
  if (length(weight)==1) 
    if (weight==-1) weight(1)=2.0; weight(2:nMB)=-1/(nMB-1); 
    else weight(1:nMB)=weight/nMB;
    end;
  end;
  for ni=1:nMB,
    n = gm.Alg.minibucket(bucket).nodes(ni);
    gm.Alg.nodes(n).wt = weight(ni);
  end;

