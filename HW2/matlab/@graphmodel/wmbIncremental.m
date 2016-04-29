function [gm lnZ] = wmbIncremental(gm, order, iBound, sBound, task, nPer)
% [gm lnZ] = wmbIncremental(gm, order, iBound, sBound, task, nPer)
% Incrementally build a weighted mini-bucket, passing nPer back/fwd messages before each merge
% TODO: fix match patterns after merges to be equivalent to inital build

[gm lnZ] = wmbInit(gm, order, 1, inf, task);   % build initial WMB bound with no cliques
% (TODO): add way to extract lnZ from wmbInit
%[gm lnZ] = wmbFwd(gm, false, false);     % unnecessary forward pass, but no matching => model unchanged
fprintf('Initial: %f\n',lnZ); 
done = false;

while (~done),
  % pass messages (nPer) to update bound at current cliques
  for it=1:nPer,
    gm=wmbBwd(gm); 
    [gm lnZnew]=wmbFwd(gm, true,false);          % update forward (param/weight updates?)
    if (strcmp(task,'sum-') && lnZ < lnZnew) lnZ=lnZnew;
    elseif (lnZ > lnZnew) lnZ=lnZnew;           % save the best bound so far
    end;
    fprintf('  it %d: %f  => %f\n',it,lnZnew,lnZ);
  end;
  [s,C] = wmbTestMerge(gm, iBound, inf);
  if (isempty(s)) 
    %wmbPrint(gm), pause;
    done = true;
  else
    [mx,idx] = max(abs(s));                    % find the largest change in the bound
%wmbPrint(gm); pause;
    fprintf('Merging => {'); for v=C{idx}, fprintf('%d ',v); end; fprintf('}\n');
    gm=wmbAddClique(gm,C{idx});                % add the new clique to the minibucket
    [gm,lnZnew]=wmbFwd(gm, false,false);       % no matching => should always improve (if reparam done correctly)
    fprintf('Updated: %f  (%f)\n',lnZnew,s(idx));
%wmbPrint(gm), pause;
    if (strcmp(task,'sum-') && lnZ < lnZnew) lnZ=lnZnew;
    elseif (lnZ > lnZnew) lnZ=lnZnew;
    end;
  end;
end;
fprintf('===> %f\n',lnZ);

