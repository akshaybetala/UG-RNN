function [lnZUp, lnZLo] = wmbSearch(mbUp,mbLo)



lnZUp = 0.0; 
lnZLo = 0.0;
for n=1:length(mbUp.Alg.nodes), if (mbUp.Alg.nodes(n).parent==0), lnZUp=lnZUp+table(mbUp.Alg.nodes(n).msgFwd); end; end;
for n=1:length(mbLo.Alg.nodes), if (mbLo.Alg.nodes(n).parent==0), lnZLo=lnZLo+table(mbLo.Alg.nodes(n).msgFwd); end; end;

% lnZ priority:  log( exp(up)-exp(lo) ) = log( exp(up) * (1-exp(lo)/exp(up)) ) = 
calcPriority = @(up,lo) up + log(1-exp(min(lo,up)-up));

order = mbUp.Alg.order;
nV    = length(order);
dims  = mbUp.dims;

nodes(1).x=[];
nodes(1).v=[];
nodes(1).hUp=lnZUp;
nodes(1).hLo=lnZLo;
nodes(1).priority = calcPriority(lnZUp,lnZLo);

freeList=[]; nFree=0;

UpIt=[]; LoIt=[];
for expanded=1:100000,
  [bestp best] = max( [nodes(:).priority] );
  if (bestp == -inf) fprintf('\nSearch Complete\n'); break; end;
  depth=length(nodes(best).x); 
  hUp = nodes(best).hUp; hLo = nodes(best).hLo;
  fprintf('(%f : %f) ',lnZUp,lnZLo);
    fprintf('  Expanding %d [%d] (%f = %f / %f), [',best,depth,bestp,hUp,hLo);
    fprintf('%d ',nodes(best).x); fprintf(']=['); fprintf('%d ',nodes(best).v); fprintf(']\n');
  if (depth==nV),    % numerical issues in upper/lower bounds? just remove node and continue
    nodes(best).priority = -inf; nFree=nFree+1; freeList(nFree)=best; continue;
  end;
  xi = order(nV-depth); 
  [xalpha perm] = sort([nodes(best).x xi]);
  sumUp = -inf; sumLo = -inf;
  for vi=1:dims(xi),
    hUpi=hUp; hLoi=hLo;
    if (nFree>0), n=freeList(nFree); nFree=nFree-1; else n=length(nodes)+1; end;
    % add node with xalpha = valpha:
    valpha = [nodes(best).v vi]; valpha = valpha(perm);
    nodes(n).x=xalpha; nodes(n).v=valpha;
%%% TODO:BUG: subtracts msgFwd too many times!!!
    for j=mbUp.Alg.minibucket(nV-depth).nodes,  % update upper bound at this node
      hUpi = hUpi + value2(mbUp.Alg.nodes(j).theta, xalpha, valpha) - value2(mbUp.Alg.nodes(j).msgFwd,xalpha,valpha);
      %hUpi = hUpi + value2(mbUp.Alg.nodes(j).theta, xalpha, valpha);
      %if (vi==1) hUpi = hUpi - value2(mbUp.Alg.nodes(j).msgFwd,xalpha,valpha); end;
      for c=mbUp.Alg.nodes(j).children, hUpi=hUpi+value2(mbUp.Alg.nodes(c).msgFwd,xalpha,valpha); end;
    end;
    for j=mbLo.Alg.minibucket(nV-depth).nodes,  % update lower bound at this node
%%% TODO:BUG: subtracts msgFwd too many times!!!
      hLoi = hLoi + value2(mbLo.Alg.nodes(j).theta, xalpha, valpha) - value2(mbLo.Alg.nodes(j).msgFwd,xalpha,valpha);
      %hLoi = hLoi + value2(mbLo.Alg.nodes(j).theta, xalpha, valpha);
      %if (vi==1) hLoi = hLoi - value2(mbLo.Alg.nodes(j).msgFwd,xalpha,valpha); end;
      for c=mbLo.Alg.nodes(j).children, hLoi=hLoi+value2(mbLo.Alg.nodes(c).msgFwd,xalpha,valpha); end;
    end;
    nodes(n).x=xalpha; nodes(n).v=valpha; nodes(n).hUp=hUpi; nodes(n).hLo=hLoi; nodes(n).priority=calcPriority(hUpi,hLoi);
    if (vi==1), sumUp=hUpi; else sumUp = sumUp + log(1 + exp(hUpi-sumUp)); end;
    if (vi==1), sumLo=hLoi; else sumLo = sumLo + log(1 + exp(hLoi-sumLo)); end;
  end;
  lnZUp = lnZUp + log(1 - exp(hUp-lnZUp) + exp(sumUp-lnZUp)); %=log( exp(lnZUp) - exp(hUp) + exp(sumUp) );
  lnZLo = lnZLo + log(1 - exp(hLo-lnZLo) + exp(sumLo-lnZLo)); %=log( exp(lnZUp) - exp(hUp) + exp(sumUp) );
  nodes(best).priority = -inf; nFree=nFree+1; freeList(nFree)=best;


  if (mod(expanded,1000)==0) 
     UpIt(end+1) = lnZUp; LoIt(end+1)=lnZLo; x = 1:length(LoIt);
    figure(1); plot(x,UpIt,'b-',x,LoIt,'r-');
  end;
end; 
    




