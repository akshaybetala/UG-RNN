function samples = gibbs(gm, nSamples,nBurn,nStep,st)
% samples = gibbs(gm, nSamples [,nBurn, nStep, init]) -- perform Gibbs sampling on graphical model gm
%  nSamples : total number of samples to return
%  nBurn    : number of burn-in steps in the Markov chain (default 0)
%  nStep    : number of steps between retained samples (default 1) 
%  init     : initial state vector (size 1 x gm.nVar)

if (nargin < 5) st = ceil(rand(1,gm.nVar).*gm.dims); end;   % random initial state value
if (nargin < 4) nStep = 1; end;
if (nargin < 3) nBurn = 0; end;
if (nargin < 2) error('Need at least two inputs'); end;
samples=zeros(nSamples,gm.nVar);

%fprintf('Running gibbs sampling for %d iteratons...\n',nBurn+nSamples*nStep);

for it=1:nBurn+nSamples*nStep,
  for x=1:gm.nVar,
    if (gm.dims(x)==0) continue; end;           
    fids = withVariable(gm,x);                         % collect all functions depending on x
    F = factor(x, ones(gm.dims(x),1));                 %
    for j=1:length(fids),                              %
      fj = getFactor(gm,j); vj = vdiff(vars(fj),x);    %   and condition on x's neighbors values
      fj = condition(fj, vj, st(1,vj));                %
      F = F * fj;
    end;
    if (sum(F)==0) st(1,x)=ceil(rand(1)*gm.dims(x));   % invalid state => move randomly
    else st(1,x) = sample(F,1);                        % otherwise draw gibbs sample
    end;
  end;
  s=fix((it-nBurn)/nStep);      % compute sample number
  if (s>0 && it==s*nStep+nBurn) % if past burn-in and last in step
    samples(s,:) = st;          %   save the state of the markov chain
  end;
end;

      

