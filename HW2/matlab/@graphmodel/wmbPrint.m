function wmbPrint(gm)
% wmbPrint(gm) : print WMB data structure (clique & joingraph) info for WMB algorithm

if (~strcmp(gm.Alg.name,'WMB'))
  fprintf('Graphical model has not been specialized to the WMB algorithm\n');
  return;
end;

for i=1:length(gm.Alg.minibucket),
  fprintf('%03d: ',gm.Alg.minibucket(i).var);
  for c=gm.Alg.minibucket(i).nodes,
    fprintf('%d {',c); 
    for v=gm.Alg.nodes(c).clique, fprintf('%d ',v); end;
    fprintf('},%0.2f => %d ; ',gm.Alg.nodes(c).wt, gm.Alg.nodes(c).parent);
  end;
  fprintf('\n');
end;


